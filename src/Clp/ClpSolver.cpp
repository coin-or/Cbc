// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Need these up front to define symbols for other imports
#include "ClpConfig.h"
#include "CoinUtilsConfig.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#define GLP_UNDEF 1
#define GLP_FEAS 2
#define GLP_INFEAS 3
#define GLP_NOFEAS 4
#define GLP_OPT 5

#if PRICE_USE_OPENMP
#include "omp.h"
#endif
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#include "CoinHelperFunctions.hpp"
#include "CoinParam.hpp"
#include "CoinPragma.hpp"
#include "CoinSort.hpp"
// History since 1.0 at end
#include "CoinFileIO.hpp"
#include "CoinModel.hpp"
#include "CoinMpsIO.hpp"
#include "CoinSignal.hpp"
#include "CoinTime.hpp"
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
#include "ClpParam.hpp"
#include "ClpParamUtils.hpp"
#include "ClpParameters.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpPresolve.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
#include "ClpSolver.hpp"
#include "ClpOutput.hpp"

#include "ClpModelParameters.hpp"

//#############################################################################
//#############################################################################

void printGeneralMessage(ClpSimplex &model, std::string message, int type)
{
  if (message.length()) {
    model.messageHandler()->message(type, model.messages())
        << message << CoinMessageEol;
  }
}

//#############################################################################
//#############################################################################

void printGeneralWarning(ClpSimplex &model, std::string message, int type)
{
  if (message.length()) {
    model.messageHandler()->message(type, model.messages())
        << message << CoinMessageEol;
  }
}

//#############################################################################
//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

CLPLIB_EXPORT
void ClpMain0(ClpSimplex &model)
{
  model.setPerturbation(50);
  model.messageHandler()->setPrefix(false);
#if CLP_INHERIT_MODE > 1
  model.setDualTolerance(1.0e-6);
  model.setPrimalTolerance(1.0e-6);
#endif
}

//#############################################################################
//#############################################################################
// old way
CLPLIB_EXPORT
void ClpMain0(ClpSimplex *model)
{
  model->setPerturbation(50);
  model->messageHandler()->setPrefix(false);
#if CLP_INHERIT_MODE > 1
  model->setDualTolerance(1.0e-6);
  model->setPrimalTolerance(1.0e-6);
#endif
}
CLPLIB_EXPORT
int ClpMain1(int argc, const char *argv[], ClpSimplex *model)
{
  std::deque<std::string> inputQueue;
  CoinParamUtils::formInputQueue(inputQueue, "clp", argc, const_cast< char ** >(argv));
  return ClpMain1(inputQueue,*model);
}
CLPLIB_EXPORT
int ClpMain1(std::deque<std::string> inputQueue, ClpSimplex &model,
             ampl_info *info)
{
  std::ostringstream buffer;
  std::string field, message, fileName;
  FILE *fp;
  ClpSimplex &model_ = model;

  // Set up all non-standard stuff
  // int numberModels=1;
#ifdef CLP_USEFUL_PRINTOUT
  double startElapsed = CoinGetTimeOfDay();
  double startCpu = CoinCpuTime();
  static std::string mpsFile = "";
  memset(debugInt, 0, sizeof(debugInt));
  memset(debugDouble, 0, sizeof(debugDouble));
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
  int doCrash = 0;
  int doVector = 0;
  int doSprint = -1;
  // set reasonable defaults
#if CLP_INHERIT_MODE > 1
#define DEFAULT_PRESOLVE_PASSES 20
#else
#define DEFAULT_PRESOLVE_PASSES 10
#endif
  int preSolve = DEFAULT_PRESOLVE_PASSES;
  bool preSolveFile = false;
  const char dirsep = CoinFindDirSeparator();
  int basisHasValues = 0;
  int substitution = 3;
  int dualize = 3; // dualize if looks promising
  ClpParameters parameters;
  parameters.setModel(&model);
  parameters[ClpParam::DUALBOUND]->setVal(model_.dualBound());
  parameters[ClpParam::DUALTOLERANCE]->setVal(model_.dualTolerance());
  parameters[ClpParam::IDIOT]->setVal(doIdiot);
  parameters[ClpParam::LOGLEVEL]->setVal(model_.logLevel());
  parameters[ClpParam::MAXFACTOR]->setVal(model_.factorizationFrequency());
  parameters[ClpParam::MAXITERATION]->setVal(model_.maximumIterations());
  parameters[ClpParam::OUTPUTFORMAT]->setVal(outputFormat);
  parameters[ClpParam::PRESOLVEPASS]->setVal(preSolve);
  parameters[ClpParam::PERTVALUE]->setVal(model_.perturbation());
  parameters[ClpParam::PRIMALTOLERANCE]->setVal(model_.primalTolerance());
  parameters[ClpParam::PRIMALWEIGHT]->setVal(model_.infeasibilityCost());
  parameters[ClpParam::TIMELIMIT]->setVal(model_.maximumSeconds());
  parameters[ClpParam::SPRINT]->setVal(doSprint);
  parameters[ClpParam::SUBSTITUTION]->setVal(substitution);
  parameters[ClpParam::DUALIZE]->setVal(dualize);
  parameters[ClpParam::PRESOLVETOLERANCE]->setVal(1.0e-8);
  int verbose = 0;

  // total number of commands read
  int numberGoodCommands = 0;
  bool goodModel = false;
  if (model_.numberRows() || model_.numberColumns()) {
    // model already built
    goodModel = true;
    numberGoodCommands = 1;
  }

  bool usingAmpl = false;
  CoinMessageHandler *generalMessageHandler = model_.messageHandler();
  generalMessageHandler->setPrefix(false);
  CoinMessages generalMessages = model_.messages();

  if (info) {
    // We're using AMPL
    usingAmpl = true;
    parameters[ClpParam::LOGLEVEL]->setVal(info->logLevel);
    goodModel = true;
  }

  // Hidden stuff for barrier
  int choleskyType = 0;
  int gamma = 0;
  parameters[ClpParam::BARRIERSCALE]->setVal(2);
  int scaleBarrier = 2;
  int doKKT = 0;
  int crossover = 2; // do crossover unless quadratic
  bool canOpen;

  // model_.scaling(1);
  // model_.setDualBound(1.0e6);
  // model_.setDualTolerance(1.0e-7);
  // ClpDualRowSteepest steep;
  // model_.setDualRowPivotAlgorithm(steep);
  // ClpPrimalColumnSteepest steepP;
  // model_.setPrimalColumnPivotAlgorithm(steepP);
  int status, iValue;
  double dValue;

  // If no arguments, print help and exit (unless called from AMPL)
  if (inputQueue.empty()) {
    if (!usingAmpl) {
      std::cout << "Clp version " << CLP_VERSION
                << " — COIN-OR Linear Programming solver\n\n"
                << "Usage:\n"
                << "  clp <model.mps[.gz]> [options] -solve\n\n"
                << "Use 'clp -help' for a list of parameters.\n";
      return 0;
    }
  } else {
    // See if first is file
    std::string inputFile = inputQueue.front();
    if (inputFile[0]!='-') {
      FILE * fp = fopen(inputFile.c_str(),"r");
      if (fp) {
	fclose(fp);
	// insert -import
	inputQueue.push_front("-import");
      }
    }
  }

  while (1) {

    // get next command
    field = CoinParamUtils::getNextField(inputQueue);

    // exit if null or similar
    if (!field.length()) {
      if (numberGoodCommands == 1 && goodModel) {
        // we just had file name - do dual or primal
        field = "-either";
      } else {
        break;
      }
    }

    // Handle field prefix stripping and special tokens
    if (field == "-") {
      // '-' alone is no longer valid; ignore and continue
      continue;
    } else if (field[0] != '-') {
       // special dispensation - taken as -import name, put name back on queue
       inputQueue.push_front(field);
       field = "import";
    } else {
       if (field != "--") {
          // take off -
          field = field.substr(1);
       } else {
          // special dispensation - taken as -import --
          field = "import";
       }
    }

    int numberMatches(0), numberShortMatches(0), numberQuery(0);

    // find out if valid command
    int paramCode = CoinParamUtils::lookupParam(field,
                                                parameters.paramVec(),
                                                &numberMatches,
                                                &numberShortMatches,
                                                &numberQuery);

    if (numberQuery > 0 || (numberShortMatches > 0 && !numberMatches)){
       continue;
    }
    if (!numberMatches) {
       std::cout << "Unrecognized parameter - " << field
       << ", exiting..."
       << std::endl;
       paramCode = ClpParam::EXIT;
    }

    // Do some translation for backwards compatibility
    switch (paramCode){
     case ClpParam::READMODEL_OLD:
       paramCode = ClpParam::READMODEL;
       break;
     case ClpParam::WRITEGMPLSOL_OLD:
       paramCode = ClpParam::WRITEGMPLSOL;
       break;
     case ClpParam::WRITEMODEL_OLD:
       paramCode = ClpParam::WRITEMODEL;
       break;
     case ClpParam::WRITESOL_OLD:
       paramCode = ClpParam::WRITESOL;
       break;
     case ClpParam::WRITESOLBINARY_OLD:
       paramCode = ClpParam::WRITESOLBINARY;
       break;
     default:
       break;
    }

    ClpParam *param = parameters[paramCode];

    int status;
    numberGoodCommands++;
    if (paramCode == ClpParam::GENERALQUERY ||
	paramCode == ClpParam::FULLGENERALQUERY) {
      // TODO Make this a method in the settings class
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
      int maxAcross = 10;
      // bool evenHidden = false;
      int commandPrintLevel =
          parameters[ClpParam::COMMANDPRINTLEVEL]->modeVal();
      if (paramCode == ClpParam::FULLGENERALQUERY) {
	verbose = 1;
	commandPrintLevel = 1;
      }
      if ((verbose & 8) != 0) {
        // even hidden
        // evenHidden = true;
        verbose &= ~8;
      }
      if (verbose < 4 && usingAmpl){
        verbose += 4;
      }
      if (verbose){
        maxAcross = 1;
      }
      std::vector<std::string> types;
      types.push_back("Invalid parameters:");
      types.push_back("Action parameters:");
      types.push_back("Integer parameters:");
      types.push_back("Double parameters:");
      types.push_back("String parameters:");
      types.push_back("Directory parameters:");
      types.push_back("File parameters:");
      types.push_back("Keyword parameters:");
      // correct types
      for (int type = 1; type < 8; type++) {
	int across = 0;
	int lengthLine = 0;
	bool first = true;
	for (int iParam = ClpParam::FIRSTPARAM + 1; iParam < ClpParam::LASTPARAM;
	     iParam++) {
	  ClpParam *p = parameters[iParam];
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
          int length = p->matchName().length();
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
      std::cout << std::endl;
      continue;
    }

    if (param->type() == CoinParam::paramDbl) {
      // get next field as double
      if (status = param->readValue(inputQueue, dValue, &message)){
        printGeneralMessage(model_, message);
        continue;
      }
      if (param->setVal(dValue, &message)){
         printGeneralMessage(model_, message);
         continue;
      } else {
         printGeneralMessage(model_, message, CLP_GENERAL2); // param-change info, suppress at default log level
      }
#if 0
      if (paramCode == ClpParam::DUALTOLERANCE)
         model_.setDualTolerance(dValue);
      else if (paramCode == ClpParam::PRIMALTOLERANCE)
         model_.setPrimalTolerance(dValue);
      else if (paramCode == ClpParam::ZEROTOLERANCE)
         model_.setSmallElementValue(dValue);
      else if (paramCode == ClpParam::DUALBOUND)
         model_.setDualBound(dValue);
      else if (paramCode == ClpParam::PRIMALWEIGHT)
         model_.setInfeasibilityCost(dValue);
      else if (paramCode == ClpParam::TIMELIMIT)
         model_.setMaximumSeconds(dValue);
      else if (paramCode == ClpParam::OBJSCALE)
         model_.setObjectiveScale(dValue);
      else if (paramCode == ClpParam::RHSSCALE)
         model_.setRhsScale(dValue);
      else if (paramCode == ClpParam::PRESOLVETOLERANCE)
         model_.setDblParam(ClpPresolveTolerance, dValue);
      else if (paramCode == ClpParam::PROGRESS)
         model_.setMinIntervalProgressUpdate(dValue);
#endif
    } else if (param->type() == CoinParam::paramInt) {
      // get next field as int
       if (status = param->readValue(inputQueue, iValue, &message)){
        printGeneralMessage(model_, message);
        continue;
      }
      if (param->setVal(iValue, &message)){
         printGeneralMessage(model_, message);
         continue;
      } else {
         printGeneralMessage(model_, message, CLP_GENERAL2); // param-change info, suppress at default log level
      }
      if (paramCode == ClpParam::PRESOLVEPASS)
         preSolve = iValue;
      else if (paramCode == ClpParam::IDIOT)
         doIdiot = iValue;
      else if (paramCode == ClpParam::SPRINT)
         doSprint = iValue;
      else if (paramCode == ClpParam::OUTPUTFORMAT)
         outputFormat = iValue;
      else if (paramCode == ClpParam::SLPVALUE)
         slpValue = iValue;
      else if (paramCode == ClpParam::CPP)
         cppValue = iValue;
      else if (paramCode == ClpParam::PRESOLVEOPTIONS)
         presolveOptions = iValue;
      else if (paramCode == ClpParam::PRINTOPTIONS)
         printOptions = iValue;
      else if (paramCode == ClpParam::SUBSTITUTION)
         substitution = iValue;
      else if (paramCode == ClpParam::DUALIZE)
         dualize = iValue;
      else if (paramCode == ClpParam::VERBOSE)
         verbose = iValue;
#if 0
      else if (paramCode == ClpParam::MAXFACTOR)
         model_.factorization()->maximumPivots(iValue);
      else if (paramCode == ClpParam::PERTVALUE)
         model_.setPerturbation(iValue);
      else if (paramCode == ClpParam::MAXITERATION)
         model_.setMaximumIterations(iValue);
      else if (paramCode == ClpParam::SPECIALOPTIONS)
         model_.setSpecialOptions(iValue);
      else if (paramCode == ClpParam::RANDOMSEED)
         model_.setRandomSeed(iValue);
      else if (paramCode == ClpParam::MORESPECIALOPTIONS)
         model_.setMoreSpecialOptions(iValue);
      else if (paramCode == ClpParam::VECTOR_MODE)
         model_.setVectorMode(iValue);
#endif
    } else if (param->type() == CoinParam::paramKwd) {
      // one of several strings
      if (status = param->readValue(inputQueue, field, &message)){
        printGeneralMessage(model_, message);
        continue;
      }
      if (param->setVal(field, &message)) {
         printGeneralMessage(model_, message);
         continue;
      }
      // Note: successful setVal for keyword params currently gives no message, but suppress if it does
      int mode = param->modeVal();
      // TODO this should be part of the push method
      switch (paramCode) {
        case ClpParam::DIRECTION:
          if (mode == 0) {
            model_.setOptimizationDirection(1);
          } else if (mode == 1) {
            model_.setOptimizationDirection(-1);
          } else {
            model_.setOptimizationDirection(0);
          }
          break;
        case ClpParam::DUALPIVOT:
          if (mode == 0) {
            ClpDualRowSteepest steep(3);
            model_.setDualRowPivotAlgorithm(steep);
          } else if (mode == 1) {
            // ClpDualRowDantzig dantzig;
            ClpDualRowDantzig dantzig;
            model_.setDualRowPivotAlgorithm(dantzig);
          } else if (mode == 2) {
            // partial steep
            ClpDualRowSteepest steep(2);
            model_.setDualRowPivotAlgorithm(steep);
          } else if (mode == 3) {
            ClpDualRowSteepest steep;
            model_.setDualRowPivotAlgorithm(steep);
          } else if (mode == 4) {
            // Positive edge steepest
            ClpPEDualRowSteepest p(fabs(parameters[ClpParam::PSI]->dblVal()));
            model_.setDualRowPivotAlgorithm(p);
          } else if (mode == 5) {
            // Positive edge Dantzig
            ClpPEDualRowDantzig p(fabs(parameters[ClpParam::PSI]->dblVal()));
            model_.setDualRowPivotAlgorithm(p);
          }
          break;
        case ClpParam::PRIMALPIVOT:
          if (mode == 0) {
            ClpPrimalColumnSteepest steep(3);
            model_.setPrimalColumnPivotAlgorithm(steep);
          } else if (mode == 1) {
            ClpPrimalColumnSteepest steep(0);
            model_.setPrimalColumnPivotAlgorithm(steep);
          } else if (mode == 2) {
            ClpPrimalColumnDantzig dantzig;
            model_.setPrimalColumnPivotAlgorithm(dantzig);
          } else if (mode == 3) {
            ClpPrimalColumnSteepest steep(4);
            model_.setPrimalColumnPivotAlgorithm(steep);
          } else if (mode == 4) {
            ClpPrimalColumnSteepest steep(1);
            model_.setPrimalColumnPivotAlgorithm(steep);
          } else if (mode == 5) {
            ClpPrimalColumnSteepest steep(2);
            model_.setPrimalColumnPivotAlgorithm(steep);
          } else if (mode == 6) {
            ClpPrimalColumnSteepest steep(10);
            model_.setPrimalColumnPivotAlgorithm(steep);
          } else if (mode == 7) {
            // Positive edge steepest
            ClpPEPrimalColumnSteepest p(
                fabs(parameters[ClpParam::PSI]->dblVal()));
            model_.setPrimalColumnPivotAlgorithm(p);
          } else if (mode == 8) {
            // Positive edge Dantzig
            ClpPEPrimalColumnDantzig p(
                fabs(parameters[ClpParam::PSI]->dblVal()));
            model_.setPrimalColumnPivotAlgorithm(p);
          }
          break;
        case ClpParam::SCALING:
          model_.scaling(mode);
          break;
        case ClpParam::AUTOSCALE:
          model_.setAutomaticScaling(mode != 0);
          break;
        case ClpParam::SPARSEFACTOR:
          model_.setSparseFactorization((1 - mode) != 0);
          break;
        case ClpParam::BIASLU:
          model_.factorization()->setBiasLU(mode);
          break;
        case ClpParam::PERTURBATION:
          if (mode == 1) {
	    if (model_.perturbation()==100)
	      model_.setPerturbation(50);
          } else {
            model_.setPerturbation(100);
	    parameters[ClpParam::PERTVALUE]->setVal(100);
	  }
          break;
        case ClpParam::ERRORSALLOWED:
          allowImportErrors = mode;
          break;
        case ClpParam::ABCWANTED:
#if   PRICE_USE_OPENMP
          omp_set_num_threads(mode);
#endif
          break;
        case ClpParam::INTPRINT:
          printMode = mode;
          break;
        case ClpParam::KEEPNAMES:
          keepImportNames = mode;
          break;
        case ClpParam::PRESOLVE:
          if (mode == 0)
            preSolve = DEFAULT_PRESOLVE_PASSES;
          else if (mode == 1)
            preSolve = 0;
          else if (mode == 2)
            preSolve = 10;
          else
            preSolveFile = true;
          break;
        case ClpParam::PFI:
          model_.factorization()->setForrestTomlin(mode == 0);
          break;
        case ClpParam::FACTORIZATION:
          model_.factorization()->forceOtherFactorization(mode);
          break;
        case ClpParam::CRASH:
          doCrash = mode;
          break;
        case ClpParam::VECTOR:
          doVector = mode;
          break;
        case ClpParam::MESSAGES:
          model_.messageHandler()->setPrefix(mode != 0);
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
    } else if (param->type() == CoinParam::paramDir){
       if (status = param->readValue(inputQueue, field, &message)){
          printGeneralMessage(model_, message);
          continue;
       }
       if (param->setVal(field, &message)){
          printGeneralMessage(model_, message);
          continue;
       }
    } else if (param->type() == CoinParam::paramFile){
       if (status = param->readValue(inputQueue, field, &message)){
          printGeneralMessage(model_, message);
          continue;
       }
       if (param->setVal(field, &message)){
          printGeneralMessage(model_, message);
          continue;
       }
    } else {
      // action
      if (paramCode == ClpParam::EXIT ||
	  paramCode == ClpParam::STOP ||
	  paramCode == ClpParam::QUIT ||
	  paramCode == ClpParam::END) {
        if (usingAmpl) {
          writeAmpl(info);
          freeArrays2(info);
          freeArgs(info);
        }
        break; // stop all
      }
#ifdef CBC_CLUMSY_CODING
      /* Synchronize Clp model - Int and Dbl */
      parameters.synchronizeModel();
#endif

      switch (paramCode) {
      case ClpParam::DUALSIMPLEX:
      case ClpParam::PRIMALSIMPLEX:
      case ClpParam::EITHERSIMPLEX:
      case ClpParam::SOLVE:
      case ClpParam::BARRIER:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          // openblas_set_num_threads(4);
          // deal with positive edge
          double psi = parameters[ClpParam::PSI]->dblVal();
          if (psi > 0.0) {
            ClpDualRowPivot *dualp = model_.dualRowPivot();
            ClpDualRowSteepest *d1 = dynamic_cast<ClpDualRowSteepest *>(dualp);
            ClpDualRowDantzig *d2 = dynamic_cast<ClpDualRowDantzig *>(dualp);
            if (d1) {
              ClpPEDualRowSteepest p(psi, d1->mode());
              model_.setDualRowPivotAlgorithm(p);
            } else if (d2) {
              ClpPEDualRowDantzig p(psi);
              model_.setDualRowPivotAlgorithm(p);
            }
            ClpPrimalColumnPivot *primalp = model_.primalColumnPivot();
            ClpPrimalColumnSteepest *p1 =
                dynamic_cast<ClpPrimalColumnSteepest *>(primalp);
            ClpPrimalColumnDantzig *p2 =
                dynamic_cast<ClpPrimalColumnDantzig *>(primalp);
            if (p1) {
              ClpPEPrimalColumnSteepest p(psi, p1->mode());
              model_.setPrimalColumnPivotAlgorithm(p);
            } else if (p2) {
              ClpPEPrimalColumnDantzig p(psi);
              model_.setPrimalColumnPivotAlgorithm(p);
            }
          }
          if (paramCode == ClpParam::EITHERSIMPLEX ||
              paramCode == ClpParam::SOLVE) {
            model_.setMoreSpecialOptions(
                16384 | model_.moreSpecialOptions());
            paramCode = ClpParam::EITHERSIMPLEX;
          }
          double objScale = parameters[ClpParam::OBJSCALE2]->dblVal();
          if (objScale != 1.0) {
            int iColumn;
            int numberColumns = model_.numberColumns();
            double *dualColumnSolution = model_.dualColumnSolution();
            ClpObjective *obj = model_.objectiveAsObject();
            assert(dynamic_cast<ClpLinearObjective *>(obj));
            double offset;
            double *objective = obj->gradient(NULL, NULL, offset, true);
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
              dualColumnSolution[iColumn] *= objScale;
              objective[iColumn] *= objScale;
              ;
            }
            int iRow;
            int numberRows = model_.numberRows();
            double *dualRowSolution = model_.dualRowSolution();
            for (iRow = 0; iRow < numberRows; iRow++)
              dualRowSolution[iRow] *= objScale;
            model_.setObjectiveOffset(objScale *
                                              model_.objectiveOffset());
          }
          ClpSolve::SolveType method;
          ClpSolve::PresolveType presolveType;
          ClpSolve solveOptions;
          ClpSimplex *model2 = &model_;
          if (paramCode == ClpParam::EITHERSIMPLEX)
            solveOptions.setSpecialOption(3, 0); // allow +-1
          if (dualize == 4) {
            solveOptions.setSpecialOption(4, 77);
            dualize = 0;
          }
          if (dualize) {
            bool tryIt = true;
            double fractionColumn = 1.0;
            double fractionRow = 1.0;
            if (dualize == 3) {
              dualize = 1;
              int numberColumns = model2->numberColumns();
              int numberRows = model2->numberRows();
              if (numberRows < 50000 || 5 * numberColumns > numberRows) {
                tryIt = false;
              } else {
                fractionColumn = 0.1;
                fractionRow = 0.3;
              }
            }
            if (tryIt) {
              ClpSimplex *thisModel =
                  static_cast<ClpSimplexOther *>(model2)->dualOfModel(
                      fractionRow, fractionColumn);
              if (thisModel) {
                buffer.str("");
                buffer << "Dual of model has " << model_.numberRows()
                       << " rows and " << model_.numberColumns()
                       << " columns" << std::endl;
                printGeneralMessage(model_, buffer.str());
                model_.setOptimizationDirection(1.0);
                model2 = thisModel;
              } else {
                thisModel = &model;
                dualize = 0;
              }
            } else {
              dualize = 0;
            }
          }
          if (preSolveFile)
            presolveOptions |= 0x40000000;
          // allow dependency
          presolveOptions |= 32768;
          solveOptions.setPresolveActions(presolveOptions);
          solveOptions.setSubstitution(substitution);
          if (preSolve != DEFAULT_PRESOLVE_PASSES && preSolve) {
            presolveType = ClpSolve::presolveNumber;
            if (preSolve < 0) {
              preSolve = -preSolve;
              if (preSolve <= 100) {
                presolveType = ClpSolve::presolveNumber;
                solveOptions.setDoSingletonColumn(true);
              } else {
                preSolve -= 100;
                presolveType = ClpSolve::presolveNumberCost;
              }
              buffer.str("");
              buffer << "Doing " << preSolve
                     << " presolve passes - picking up non-costed slacks"
                     << std::endl;
              printGeneralMessage(model_, buffer.str());
            }
          } else if (preSolve) {
            presolveType = ClpSolve::presolveOn;
          } else {
            presolveType = ClpSolve::presolveOff;
          }
          solveOptions.setPresolveType(presolveType, preSolve);
          if (paramCode == ClpParam::DUALSIMPLEX) {
            method = ClpSolve::useDual;
          } else if (paramCode == ClpParam::PRIMALSIMPLEX) {
            method = ClpSolve::usePrimalorSprint;
          } else if (paramCode == ClpParam::EITHERSIMPLEX) {
            method = ClpSolve::automatic;
            if (doCrash > 6) {
              solveOptions.setSpecialOption(6, 1, doCrash - 6);
              doCrash = 0;
            }
            if (doIdiot > 0)
              solveOptions.setSpecialOption(1, 2, doIdiot);
          } else {
            method = ClpSolve::useBarrier;
            if (doIdiot > 0)
              solveOptions.setSpecialOption(1, 2, doIdiot); // dense threshold
            if (crossover == 1) {
              method = ClpSolve::useBarrierNoCross;
            } else if (crossover == 2) {
              ClpObjective *obj = model_.objectiveAsObject();
              if (obj->type() > 1) {
                method = ClpSolve::useBarrierNoCross;
                presolveType = ClpSolve::presolveOff;
                solveOptions.setPresolveType(presolveType, preSolve);
              }
            }
          }
          solveOptions.setSolveType(method);
          solveOptions.setSpecialOption(5, printOptions & 1);
          if (doVector) {
            model_.setVectorMode(doVector);
            ClpMatrixBase *matrix = model_.clpMatrix();
            if (dynamic_cast<ClpPackedMatrix *>(matrix)) {
              ClpPackedMatrix *clpMatrix =
                  dynamic_cast<ClpPackedMatrix *>(matrix);
              clpMatrix->makeSpecialColumnCopy();
            }
          }
          if (method == ClpSolve::useDual) {
            // dual
            if (doCrash && doCrash < 7)
              solveOptions.setSpecialOption(0, 1, doCrash); // crash
            else if (doIdiot)
              solveOptions.setSpecialOption(0, 2, doIdiot); // possible idiot
          } else if (method == ClpSolve::usePrimalorSprint) {
            // primal
            // if slp turn everything off
            if (slpValue > 0) {
              doCrash = 0;
              doSprint = 0;
              doIdiot = -1;
              solveOptions.setSpecialOption(1, 10, slpValue); // slp
              method = ClpSolve::usePrimal;
            }
            if (doCrash > 6) {
              solveOptions.setSpecialOption(6, 1, doCrash - 6);
              doCrash = 0;
            }
            if (doCrash > 0) {
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
          } else if (method == ClpSolve::useBarrier ||
                     method == ClpSolve::useBarrierNoCross) {
            int barrierOptions = choleskyType;
            if (scaleBarrier) {
              if ((scaleBarrier & 1) != 0)
                barrierOptions |= 8;
              barrierOptions |= 2048 * (scaleBarrier >> 1);
            }
            if (doKKT)
              barrierOptions |= 16;
            if (gamma)
              barrierOptions |= 32 * gamma;
            if (crossover == 3)
              barrierOptions |= 256; // try presolve in crossover
            solveOptions.setSpecialOption(4, barrierOptions);
          }
          int status;
          if (cppValue >= 0) {
            // generate code
            fp = fopen("user_driver.cpp", "w");
            if (fp) {
              // generate enough to do solveOptions
              model2->generateCpp(fp);
              solveOptions.generateCpp(fp);
              fclose(fp);
              // now call generate code
              generateCode("user_driver.cpp", cppValue);
            } else {
              printGeneralMessage(model_,
                                  "Unable to open file user_driver.cpp\n");
            }
          }
#ifdef CLP_MULTIPLE_FACTORIZATIONS
          int denseCode = parameters[ClpParam::DENSE]->intVal();
          if (denseCode != -1)
            model2->factorization()->setGoDenseThreshold(denseCode);
          int smallCode = parameters[ClpParam::SMALLFACT]->intVal();
          if (smallCode != -1)
            model2->factorization()->setGoSmallThreshold(smallCode);
          model2->factorization()->goDenseOrSmall(model2->numberRows());
#endif
          // Install unified LP+Idiot+Sprint progress handlers.
          // ClpLpMsgHandler intercepts Idiot/Sprint messages (ext 30/34).
          // ClpLpEventHandler intercepts endOfIteration for LP rows.
          // Both share ClpLpPhaseState so all rows appear in one table.
          // We do NOT setLogLevel(0) — Idiot checks the model log level
          // before constructing messages; setting it to 0 silences Idiot.
          // Noisy simplex messages are suppressed by raising their detail
          // level to 2 in ClpMessage.cpp; ClpLpMsgHandler suppresses the rest.
          const int lpIterFreq = parameters[ClpParam::PROGRESSITER]->intVal();
          const double lpTimeFreq = model2->getMinIntervalProgressUpdate();
          auto lpState = std::make_shared<ClpLpPhaseState>();
          lpState->fp       = model_.messageHandler()->filePointer();
          lpState->utf8     = ClpOutput::useUtf8();
          lpState->compact  = ClpOutput::useCompact();
          lpState->logLevel = model2->logLevel();
          lpState->iterFreq = lpIterFreq;
          lpState->timeFreq = lpTimeFreq;
          lpState->origRows = model2->numberRows();
          lpState->origCols = model2->numberColumns();
          lpState->startTime = CoinWallclockTime();
          lpState->lastPrintTime = lpState->startTime;
          lpState->title = "LP solve";
          ClpLpMsgHandler   lpMsgH(lpState);
          ClpLpEventHandler lpEvtH(lpState);
          bool lpMsgOldDefault;
          CoinMessageHandler *lpSavedMsg =
            model2->pushMessageHandler(&lpMsgH, lpMsgOldDefault);
          model2->passInEventHandler(&lpEvtH);
          ClpLpEventHandler *lpProg =
            dynamic_cast<ClpLpEventHandler *>(model2->eventHandler());

          try {
            status = model2->initialSolve(solveOptions);
            // Print final LP status and close the table
            if (lpProg)
              lpProg->printFinalStatus();
            // Restore handlers
            model2->popMessageHandler(lpSavedMsg, lpMsgOldDefault);
            ClpEventHandler defaultHandler;
            model2->passInEventHandler(&defaultHandler);
            lpProg = nullptr;
            if (usingAmpl) {
              double value = model2->getObjValue() * model2->getObjSense();
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
              info->problemStatus = iStat;
              info->objValue = value;
              pos +=
                  sprintf(buf + pos, " objective %.*g", ampl_obj_prec(), value);
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
              info->columnStatus =
                  reinterpret_cast<int *>(malloc(numberColumns * sizeof(int)));
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
#ifndef NDEBUG
            // if infeasible check ray
            if (model2->status() == 1) {
              ClpSimplex *simplex = model2;
              if (simplex->ray()) {
                // make sure we use non-scaled versions
                ClpPackedMatrix *saveMatrix = simplex->swapScaledMatrix(NULL);
                double *saveScale = simplex->swapRowScale(NULL);
                // could use existing arrays
                int numberRows = simplex->numberRows();
                int numberColumns = simplex->numberColumns();
                double *farkas = new double[2 * numberColumns + numberRows];
                double *bound = farkas + numberColumns;
                double *effectiveRhs = bound + numberColumns;
                // get ray as user would
                double *ray = simplex->infeasibilityRay();
                // get farkas row
                memset(farkas, 0,
                       (2 * numberColumns + numberRows) * sizeof(double));
                simplex->transposeTimes(-1.0, ray, farkas);
                // Put nonzero bounds in bound
                const double *columnLower = simplex->columnLower();
                const double *columnUpper = simplex->columnUpper();
                int numberBad = 0;
                for (int i = 0; i < numberColumns; i++) {
                  double value = farkas[i];
                  double boundValue = 0.0;
                  if (simplex->getStatus(i) == ClpSimplex::basic) {
                    // treat as zero if small
                    if (fabs(value) < 1.0e-8) {
                      value = 0.0;
                      farkas[i] = 0.0;
                    }
                    if (value) {
                      // printf("basic %d direction %d farkas %g\n",
                      //	   i,simplex->directionOut(),value);
                      if (value < 0.0)
                        boundValue = std::max(columnLower[i], -1.0e20);
                      else
                        boundValue = std::min(columnUpper[i], 1.0e20);
                    }
                  } else if (fabs(value) > 1.0e-10) {
                    if (value < 0.0)
                      boundValue = columnLower[i];
                    else
                      boundValue = columnUpper[i];
                  }
                  bound[i] = boundValue;
                  if (fabs(boundValue) > 1.0e10)
                    numberBad++;
                }
                const double *rowLower = simplex->rowLower();
                const double *rowUpper = simplex->rowUpper();
                bool printBad = simplex->logLevel() > 3;
                // int pivotRow = simplex->spareIntArray_[3];
                // bool badPivot=pivotRow<0;
                for (int i = 0; i < numberRows; i++) {
                  double value = ray[i];
                  double rhsValue = 0.0;
                  if (simplex->getRowStatus(i) == ClpSimplex::basic) {
                    // treat as zero if small
                    if (fabs(value) < 1.0e-8) {
                      value = 0.0;
                      ray[i] = 0.0;
                    }
                    if (value) {
                      // printf("row basic %d direction %d ray %g\n",
                      //	   i,simplex->directionOut(),value);
                      if (value < 0.0)
                        rhsValue = rowLower[i];
                      else
                        rhsValue = rowUpper[i];
                    }
                  } else if (fabs(value) > 1.0e-10) {
                    if (value < 0.0)
                      rhsValue = rowLower[i];
                    else
                      rhsValue = rowUpper[i];
                  }
                  effectiveRhs[i] = rhsValue;
                  if (fabs(effectiveRhs[i]) > 1.0e10 && printBad) {
                    buffer.str("");
                    buffer << "Large rhs row " << i << " " << effectiveRhs[i]
                           << " after" << std::endl;
                    printGeneralWarning(model_, buffer.str());
                  }
                }
                simplex->times(-1.0, bound, effectiveRhs);
                double bSum = 0.0;
                for (int i = 0; i < numberRows; i++) {
                  bSum += effectiveRhs[i] * ray[i];
                  if (fabs(effectiveRhs[i]) > 1.0e10 && printBad) {
                    buffer.str("");
                    buffer << "Large rhs row " << i << " " << effectiveRhs[i]
                           << " after" << std::endl;
                    printGeneralWarning(model_, buffer.str());
                  }
                }
                if ((numberBad || bSum > 1.0e-6) && printBad) {
                  buffer.str("");
                  buffer << "Bad infeasibility ray " << bSum << "  - "
                         << numberBad << " bad" << std::endl;
                  printGeneralWarning(model_, buffer.str());
                } else {
                  // printf("Good ray - infeasibility %g\n",
                  //     -bSum);
                }
                delete[] ray;
                delete[] farkas;
                simplex->swapRowScale(saveScale);
                simplex->swapScaledMatrix(saveMatrix);
              } else {
                // printf("No dual ray\n");
              }
            }
#endif
          } catch (CoinError e) {
            // Clean up progress handler on error
            if (lpProg) {
              ClpEventHandler defaultHandler;
              model2->passInEventHandler(&defaultHandler);
              lpProg = nullptr;
            }
            e.print();
            status = -1;
          }
          if (dualize) {
            ClpSimplex *thisModel = &model_;
            int returnCode =
                static_cast<ClpSimplexOther *>(thisModel)->restoreFromDual(
                    model2);
            if (model2->status() == 3)
              returnCode = 0;
            delete model2;
            if (returnCode && dualize != 2) {
              currentModel = &model_;
              // register signal handler
              signal(SIGINT, signal_handler);
              model_.primal(1);
              currentModel = NULL;
            }
            buffer.str("");
            buffer
                << "After translating dual back to primal - objective value is "
                << model_.objectiveValue();
            printGeneralMessage(model_, buffer.str());
            // switch off (user can switch back on)
            parameters[ClpParam::DUALIZE]->setVal(dualize);
          }
          if (status >= 0)
            basisHasValues = 1;
        } break;
      case ClpParam::STATISTICS:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          // If presolve on look at presolved
          bool deleteModel2 = false;
          ClpSimplex *model2 = &model_;
          if (preSolve) {
            ClpPresolve pinfo;
            int presolveOptions2 = presolveOptions & ~0x40000000;
            if ((presolveOptions2 & 0xffff) != 0)
              pinfo.setPresolveActions(presolveOptions2);
            pinfo.setSubstitution(substitution);
            if ((printOptions & 1) != 0)
              pinfo.statistics();
            double presolveTolerance =
                parameters[ClpParam::PRESOLVETOLERANCE]->dblVal();
            model2 = pinfo.presolvedModel(model_, presolveTolerance,
                                          true, preSolve);
            if (model2) {
              printGeneralMessage(model_, "Statistics for presolved model\n");
              deleteModel2 = true;
            } else {
              printGeneralMessage(model_,
                                  "Presolved model looks infeasible - will use "
                                  "unpresolved\n");
              model2 = &model_;
            }
          } else {
            printGeneralMessage(model_, "Statistics for unpresolved model\n");
            model2 = &model_;
          }
          statistics(&model_, model2);
          if (deleteModel2)
            delete model2;
        }
        break;
      case ClpParam::TIGHTEN:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          int numberInfeasibilities = model_.tightenPrimalBounds();
          if (numberInfeasibilities)
            printGeneralWarning(model_,
                                "** Analysis indicates model infeasible\n");
        } break;
      case ClpParam::PLUSMINUS:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          ClpMatrixBase *saveMatrix = model_.clpMatrix();
          ClpPackedMatrix *clpMatrix =
              dynamic_cast<ClpPackedMatrix *>(saveMatrix);
          if (clpMatrix) {
            ClpPlusMinusOneMatrix *newMatrix =
                new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
            if (newMatrix->getIndices()) {
              model_.replaceMatrix(newMatrix);
              delete saveMatrix;
              printGeneralMessage(model_,
                                  "Matrix converted to +- one matrix\n");
            } else {
              printGeneralWarning(
                  model_, "Matrix cannot be converted to +- 1 matrix\n");
            }
          } else {
            printGeneralWarning(model_, "Matrix not a ClpPackedMatrix\n");
          }
        } break;
      case ClpParam::NETWORK: {
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          ClpMatrixBase *saveMatrix = model_.clpMatrix();
          ClpPackedMatrix *clpMatrix =
              dynamic_cast<ClpPackedMatrix *>(saveMatrix);
          if (clpMatrix) {
            ClpNetworkMatrix *newMatrix =
                new ClpNetworkMatrix(*(clpMatrix->matrix()));
            if (newMatrix->getIndices()) {
              model_.replaceMatrix(newMatrix);
              delete saveMatrix;
              printGeneralMessage(model_,
                                  "Matrix converted to network matrix\n");
            } else {
              printGeneralWarning(
                  model_, "Matrix can not be converted to network matrix\n");
            }
          } else {
            printGeneralWarning(model_, "Matrix not a ClpPackedMatrix\n");
          }
        } break;
      case ClpParam::IMPORT: {
        // get next field
        status = param->readValue(inputQueue, fileName, &message);
        canOpen = false;
        bool absolutePath = true;
        CoinParamUtils::processFile(fileName,
                             parameters[ClpParam::DIRECTORY]->dirName());
        if (fileName == ""){
           fileName = parameters[ClpParam::IMPORTFILE]->fileName();
        }else{
           parameters[ClpParam::IMPORTFILE]->setFileName(fileName);
        }
        if (fileName[0] != '/' && fileName[0] != '\\' &&
            !strchr(fileName.c_str(), ':')) {
           absolutePath = false;
        }
        // See if gmpl file
        int gmpl = 0;
        std::string gmplData;
        if (field[0] == '-') {
          // stdin
          canOpen = true;
          fileName = "-";
          if (field == "-lp")
            gmpl = -1;
        } else {
          // See if .lp
          const char *c_name = fileName.c_str();
          size_t length = strlen(c_name);
          if ((length > 3 && !strncmp(c_name + length - 3, ".lp", 3)) ||
              (length > 6 && !strncmp(c_name + length - 6, ".lp.gz", 6)) ||
              (length > 7 && !strncmp(c_name + length - 7, ".lp.bz2", 7))){
             gmpl = -1; // .lp
          }
          length = fileName.size();
          size_t percent = fileName.find('%');
          if (percent < length && percent > 0) {
             printf("Clp was not built with GMPL support. Exiting.\n");
             // This is surely not the right thing to do here. Should we
             // throw an exceptioon? Exit?
             abort();
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
                buffer << "Unable to open file " << gmplData << std::endl;
                printGeneralWarning(model_, buffer.str());
                continue;
              }
            }
          } else {
            buffer.str("");
            buffer << "Unable to open file " << gmplData << std::endl;
            printGeneralWarning(model_, buffer.str());
            continue;
          }
        }
        int status;
#ifdef CLP_USEFUL_PRINTOUT
        mpsFile = fileName;
#endif
        if (!gmpl) {
           status = model_.readMps(
                fileName.c_str(), keepImportNames != 0, allowImportErrors != 0);
        } else if (gmpl > 0) {
           printGeneralWarning(
                model_, "Clp was not built with GMPL support. Exiting.\n");
            // This is surely not the right thing to do here. Should we
            // throw an exceptioon? Exit?
           abort();
        } else {
#ifdef KILL_ZERO_READLP
           status = model_.readLp(
                fileName.c_str(), model_.getSmallElementValue());
#else
           status = model_.readLp(fileName.c_str(), 1.0e-12);
#endif
        }
        if (!status || (status > 0 && allowImportErrors)) {
           goodModel = true;
           // sets to all slack (not necessary?)
           model_.createStatus();
           // Print compact problem summary (dimensions + coefficient ranges)
           ClpOutput::printProblemSummary(model_.messageHandler(), model_,
             model_.logLevel());
           // Go to canned file if just input file
           if (inputQueue.empty()) {
              // only if ends .mps
              char *find = const_cast<char *>(strstr(fileName.c_str(), ".mps"));
              if (find && find[4] == '\0') {
                 find[1] = 'p';
                 find[2] = 'a';
                 find[3] = 'r';
                 std::ifstream ifs(fileName.c_str());
                 if (ifs.is_open()) {
                    CoinParamUtils::readFromStream(inputQueue, ifs);
                 } else {
		   // NOT an error
		   //buffer.str("");
		   //buffer << "No parameter file " << fileName << " found"
		   //      << std::endl;
		   //printGeneralMessage(model_, buffer.str());
                 }
              }
           }
        } else {
           // errors
           buffer.str("");
           buffer << "There were " << status << " errors on input"
                  << std::endl;
           printGeneralMessage(model_, buffer.str());
        }
        } break;
      case ClpParam::EXPORT:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          double objScale = parameters[ClpParam::OBJSCALE2]->dblVal();
          if (objScale != 1.0) {
            int iColumn;
            int numberColumns = model_.numberColumns();
            double *dualColumnSolution = model_.dualColumnSolution();
            ClpObjective *obj = model_.objectiveAsObject();
            assert(dynamic_cast<ClpLinearObjective *>(obj));
            double offset;
            double *objective = obj->gradient(NULL, NULL, offset, true);
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
              dualColumnSolution[iColumn] *= objScale;
              objective[iColumn] *= objScale;
              ;
            }
            int iRow;
            int numberRows = model_.numberRows();
            double *dualRowSolution = model_.dualRowSolution();
            for (iRow = 0; iRow < numberRows; iRow++)
              dualRowSolution[iRow] *= objScale;
            model_.setObjectiveOffset(objScale *
                                              model_.objectiveOffset());
          }
          // get next field
          param->readValue(inputQueue, fileName, &message);
          CoinParamUtils::processFile(fileName,
                                 parameters[ClpParam::DIRECTORY]->dirName());
          if (fileName == ""){
             fileName = parameters[ClpParam::EXPORTFILE]->fileName();
          }else{
             parameters[ClpParam::EXPORTFILE]->setFileName(fileName);
          }
          fp = fopen(fileName.c_str(), "w");
          if (fp) {
            // can open - lets go for it
            fclose(fp);
          } else {
            buffer.str("");
            buffer << "Unable to open file " << fileName << std::endl;
            printGeneralWarning(model_, buffer.str());
            continue;
          }
          // If presolve on then save presolved
          bool deleteModel2 = false;
          ClpSimplex *model2 = &model_;
          if (dualize && dualize < 3) {
             model2 = static_cast<ClpSimplexOther *>(model2)->dualOfModel();
             buffer.str("");
             buffer << "Dual of model has " << model2->numberRows()
                    << " rows and " << model2->numberColumns() << " columns"
                    << std::endl;
             printGeneralMessage(model_, buffer.str());
             model2->setOptimizationDirection(1.0);
             preSolve = 0; // as picks up from model
          }
          if (preSolve) {
             ClpPresolve pinfo;
             int presolveOptions2 = presolveOptions & ~0x40000000;
             if ((presolveOptions2 & 0xfffff) != 0)
                pinfo.setPresolveActions(presolveOptions2);
             pinfo.setSubstitution(substitution);
             if ((printOptions & 1) != 0)
                pinfo.statistics();
             double presolveTolerance =
                parameters[ClpParam::PRESOLVETOLERANCE]->dblVal();
             model2 = pinfo.presolvedModel(model_, presolveTolerance,
                                           true, preSolve, false, false);
             if (model2) {
                buffer.str("");
                buffer << "Saving presolved model on " << fileName << std::endl;
                printGeneralMessage(model_, buffer.str());
                deleteModel2 = true;
             } else {
                buffer.str("");
                buffer
                   << "Presolved model looks infeasible - saving original on "
                   << fileName << std::endl;
                printGeneralMessage(model_, buffer.str());
                deleteModel2 = false;
                model2 = &model_;
             }
          } else {
             buffer.str("");
             buffer << "Saving model on " << fileName << std::endl;
             printGeneralMessage(model_, buffer.str());
          }
#if 0
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
#ifdef STRIPBLANKS
                char * xx = rowNames[iRow];
                int i;
                int length = strlen(xx);
                int n = 0;
                for (i = 0; i < length; i++) {
                   if (xx[i] != ' ')
                      xx[n++] = xx[i];
                }
                xx[n] = '\0';
#endif
             }

             columnNames = new char * [numberColumns];
             for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                columnNames[iColumn] =
                   CoinStrdup(model2->columnName(iColumn).c_str());
#ifdef STRIPBLANKS
                char * xx = columnNames[iColumn];
                int i;
                int length = strlen(xx);
                int n = 0;
                for (i = 0; i < length; i++) {
                   if (xx[i] != ' ')
                      xx[n++] = xx[i];
                }
                xx[n] = '\0';
#endif
             }
          }
          CoinMpsIO writer;
          writer.setMpsData(*model2->matrix(), COIN_DBL_MAX,
                            model2->getColLower(), model2->getColUpper(),
                            model2->getObjCoefficients(),
                            (const char*) 0 /*integrality*/,
                            model2->getRowLower(), model2->getRowUpper(),
                            columnNames, rowNames);
          // Pass in array saying if each variable integer
          writer.copyInIntegerInformation(model2->integerInformation());
          writer.setObjectiveOffset(model2->objectiveOffset());
          writer.writeMps(fileName.c_str(), 0, 1, 1);
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
#else
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
	    //OsiClpSolverInterface solver(model2);
	    model2->writeLp(fileName.c_str(),"");
	    fclose(fp);
	  }
#endif
          if (deleteModel2)
             delete model2;
        } break;
      case ClpParam::BASISIN:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          param->readValue(inputQueue, fileName, &message);
          CoinParamUtils::processFile(fileName,
                                parameters[ClpParam::DIRECTORY]->dirName(),
                                      &canOpen);
          if (!canOpen) {
             buffer.str("");
             buffer << "Unable to open file " << fileName.c_str();
             printGeneralMessage(model_, buffer.str());
             continue;
          }
          if (fileName == ""){
             fileName = parameters[ClpParam::BASISFILE]->fileName();
          }else{
             parameters[ClpParam::BASISFILE]->setFileName(fileName);
          }
          int values = model_.readBasis(fileName.c_str());
          if (values == 0)
             basisHasValues = -1;
          else
             basisHasValues = 1;
        } break;
      case ClpParam::PRINTMASK:
        if (status = param->readValue(inputQueue, field, &message)){
           printGeneralMessage(model_, message);
           continue;
        }
        if (param->setVal(field, &message)){
           printGeneralMessage(model_, message);
           continue;
        }
        break;
      case ClpParam::BASISOUT:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          param->readValue(inputQueue, fileName, &message);
          CoinParamUtils::processFile(fileName,
                                parameters[ClpParam::DIRECTORY]->dirName());
          if (fileName == ""){
             fileName = parameters[ClpParam::BASISFILE]->fileName();
          }else{
             parameters[ClpParam::BASISFILE]->setFileName(fileName);
          }
          canOpen = false;
          fp = fopen(fileName.c_str(), "w");
          if (fp) {
            // can open - lets go for it
            fclose(fp);
            canOpen = true;
          } else {
            buffer.str("");
            buffer << "Unable to open file " << fileName << std::endl;
            printGeneralWarning(model_, buffer.str());
            continue;
          }
          ClpSimplex *model2 = &model_;
          model2->writeBasis(fileName.c_str(), outputFormat > 1,
                             outputFormat - 2);
        } break;
      case ClpParam::PARAMETRICS:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          param->readValue(inputQueue, fileName, &message);
          CoinParamUtils::processFile(fileName,
                                parameters[ClpParam::DIRECTORY]->dirName(),
                                      &canOpen);
          if (!canOpen) {
             buffer.str("");
             buffer << "Unable to open file " << fileName.c_str();
             printGeneralMessage(model_, buffer.str());
             continue;
          }
          if (fileName == ""){
             fileName = parameters[ClpParam::EXPORTFILE]->fileName();
          }else{
             parameters[ClpParam::EXPORTFILE]->setFileName(fileName);
          }
          ClpSimplex *model2 = &model_;
          static_cast<ClpSimplexOther *>(model2)->parametrics(fileName.c_str());
        } break;
      case ClpParam::WRITEMODEL: {
        param->readValue(inputQueue, fileName, &message);
        CoinParamUtils::processFile(fileName,
                                    parameters[ClpParam::DIRECTORY]->dirName());
        if (fileName == ""){
           fileName = parameters[ClpParam::MODELFILE]->fileName();
        }else{
           parameters[ClpParam::MODELFILE]->setFileName(fileName);
        }
        status = param->readValue(inputQueue, field);
        canOpen = false;
        fp = fopen(fileName.c_str(), "wb");
        if (fp) {
          // can open - lets go for it
          fclose(fp);
          canOpen = true;
        } else {
          buffer.str("");
          buffer << "Unable to open file " << fileName << std::endl;
          printGeneralWarning(model_, buffer.str());
          continue;
        }
        int status;
        // If presolve on then save presolved
        bool deleteModel2 = false;
        ClpSimplex *model2 = &model_;
        if (preSolve) {
           ClpPresolve pinfo;
           double presolveTolerance =
              parameters[ClpParam::PRESOLVETOLERANCE]->dblVal();
           model2 = pinfo.presolvedModel(model_, presolveTolerance,
                                         false, preSolve);
           if (model2) {
              buffer.str("");
              buffer << "Saving presolved model on " << fileName << std::endl;
              printGeneralMessage(model_, buffer.str());
              deleteModel2 = true;
           } else {
              buffer.str("");
              buffer << "Presolved model looks infeasible - saving original on "
                     << fileName << std::endl;
              printGeneralMessage(model_, buffer.str());
              deleteModel2 = false;
              model2 = &model_;
           }
        } else {
           buffer.str("");
           buffer << "Saving model on " << fileName << std::endl;
           printGeneralMessage(model_, buffer.str());
        }
        status = model2->saveModel(fileName.c_str());
        if (deleteModel2)
           delete model2;
        if (!status) {
           goodModel = true;
        } else {
           // errors
           printGeneralWarning(model_, "There were errors on output\n");
        }
      } break;
      case ClpParam::READMODEL: {
        param->readValue(inputQueue, fileName, &message);
        CoinParamUtils::processFile(fileName,
                                    parameters[ClpParam::DIRECTORY]->dirName(),
                                    &canOpen);
        if (!canOpen) {
           buffer.str("");
           buffer << "Unable to open file " << fileName.c_str();
           printGeneralMessage(model_, buffer.str());
           continue;
        }
        if (fileName == ""){
           fileName = parameters[ClpParam::MODELFILE]->fileName();
        }else{
           parameters[ClpParam::MODELFILE]->setFileName(fileName);
        }
        int status = model_.restoreModel(fileName.c_str());
        if (!status) {
           goodModel = true;
        } else {
           // errors
           printGeneralWarning(model_, "There were errors on output\n");
        }
      } break;
      case ClpParam::MAXIMIZE:
        model_.setOptimizationDirection(-1);
        break;
      case ClpParam::MINIMIZE:
        model_.setOptimizationDirection(1);
        break;
      case ClpParam::ALLSLACK:
        model_.allSlackBasis(true);
        break;
      case ClpParam::REVERSE:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          int iColumn;
          int numberColumns = model_.numberColumns();
          double *dualColumnSolution = model_.dualColumnSolution();
          ClpObjective *obj = model_.objectiveAsObject();
          assert(dynamic_cast<ClpLinearObjective *>(obj));
          double offset;
          double *objective = obj->gradient(NULL, NULL, offset, true);
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            dualColumnSolution[iColumn] = -dualColumnSolution[iColumn];
            objective[iColumn] = -objective[iColumn];
          }
          int iRow;
          int numberRows = model_.numberRows();
          double *dualRowSolution = model_.dualRowSolution();
          for (iRow = 0; iRow < numberRows; iRow++) {
            dualRowSolution[iRow] = -dualRowSolution[iRow];
          }
          model_.setObjectiveOffset(-model_.objectiveOffset());
        } break;
      case ClpParam::NETLIB_DUAL:
      case ClpParam::NETLIB_EITHER:
      case ClpParam::NETLIB_BARRIER:
      case ClpParam::NETLIB_PRIMAL:
      case ClpParam::NETLIB_TUNE: {
        // create fields for unitTest
        const char *fields[4];
        int nFields = 4;
        fields[0] = "fake main from unitTest";
        std::string mpsfield = "-dirSample=";
        mpsfield += parameters[ClpParam::DIRSAMPLE]->dirName();
        fields[1] = mpsfield.c_str();
        std::string netfield = "-dirNetlib=";
        netfield += parameters[ClpParam::DIRNETLIB]->dirName();
        fields[2] = netfield.c_str();
        fields[3] = "-netlib";
        int algorithm;
        if (paramCode == ClpParam::NETLIB_DUAL) {
          std::cerr << "Doing netlib with dual algorithm" << std::endl;
          algorithm = 0;
          model_.setMoreSpecialOptions(
              model_.moreSpecialOptions() | 32768);
        } else if (paramCode == ClpParam::NETLIB_BARRIER) {
          std::cerr << "Doing netlib with barrier algorithm" << std::endl;
          algorithm = 2;
        } else if (paramCode == ClpParam::NETLIB_EITHER) {
          std::cerr << "Doing netlib with dual or primal algorithm"
                    << std::endl;
          algorithm = 3;
        } else if (paramCode == ClpParam::NETLIB_TUNE) {
          std::cerr << "Doing netlib with best algorithm!" << std::endl;
          algorithm = 5;
          // uncomment next to get active tuning
          // algorithm=6;
        } else {
          std::cerr << "Doing netlib with primal algorithm" << std::endl;
          algorithm = 1;
        }
        // int specialOptions = model_.specialOptions();
        model_.setSpecialOptions(0);
        ClpSolve solveOptions;
        ClpSolve::PresolveType presolveType;
        if (preSolve){
          presolveType = ClpSolve::presolveOn;
        } else {
          presolveType = ClpSolve::presolveOff;
        }
        solveOptions.setPresolveType(presolveType, 5);
        if (doSprint >= 0 || doIdiot >= 0) {
          if (doSprint > 0) {
            // sprint overrides idiot
            solveOptions.setSpecialOption(1, 3, doSprint); // sprint
          } else if (doIdiot > 0) {
            solveOptions.setSpecialOption(1, 2, doIdiot); // idiot
          } else {
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
        }
#if FACTORIZATION_STATISTICS
        {
          extern int loSizeX;
          extern int hiSizeX;
          for (int i = 0; i < argc; i++) {
            if (!strcmp(argv[i], "-losize")) {
              int size = atoi(argv[i + 1]);
              if (size > 0)
                loSizeX = size;
            }
            if (!strcmp(argv[i], "-hisize")) {
              int size = atoi(argv[i + 1]);
              if (size > loSizeX)
                hiSizeX = size;
            }
          }
          if (loSizeX != -1 || hiSizeX != 1000000) {
            buffer.str("");
            buffer << "Solving problems " << loSizeX << "<= and <" << hiSizeX
                   << std::endl;
            printGeneralMessage(model_, buffer.str());
          }
        }
#endif
          // for moment then back to model_
          int specialOptions = model_.specialOptions();
          mainTest(nFields, fields, algorithm, model_, solveOptions,
                   specialOptions, doVector != 0);
        } break;
      case ClpParam::UNITTEST: {
        // create fields for unitTest
        const char *fields[2];
        int nFields = 2;
        fields[0] = "fake main from unitTest";
        std::string dirfield = "-dirSample=";
        dirfield += parameters[ClpParam::DIRSAMPLE]->dirName();
        fields[1] = dirfield.c_str();
        int specialOptions = model_.specialOptions();
        model_.setSpecialOptions(0);
        int algorithm = -1;
        if (model_.numberRows())
          algorithm = 7;
        ClpSolve solveOptions;
        ClpSolve::PresolveType presolveType;
        if (preSolve)
          presolveType = ClpSolve::presolveOn;
        else
          presolveType = ClpSolve::presolveOff;
        solveOptions.setPresolveType(presolveType, 5);
        mainTest(nFields, fields, algorithm, model_, solveOptions,
                 specialOptions, doVector != 0);
      } break;
      case ClpParam::FAKEBOUND:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          // get bound
          if (status = param->readValue(inputQueue, dValue, &message)){
             std::cout << "Must enter value for " << param->name()
                       << std::endl;
             continue;
          }
          buffer.str("");
          buffer << "Setting " << param->name() << " to DEBUG " << dValue
                 << std::endl;
          printGeneralMessage(model_, buffer.str());
          int iRow;
          int numberRows = model_.numberRows();
          double *rowLower = model_.rowLower();
          double *rowUpper = model_.rowUpper();
          for (iRow = 0; iRow < numberRows; iRow++) {
             // leave free ones for now
             if (rowLower[iRow] > -1.0e20 || rowUpper[iRow] < 1.0e20) {
                rowLower[iRow] = std::max(rowLower[iRow], -dValue);
                rowUpper[iRow] = std::min(rowUpper[iRow], dValue);
             }
          }
          int iColumn;
          int numberColumns = model_.numberColumns();
          double *columnLower = model_.columnLower();
          double *columnUpper = model_.columnUpper();
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
             // leave free ones for now
             if (columnLower[iColumn] > -1.0e20 ||
                 columnUpper[iColumn] < 1.0e20) {
                columnLower[iColumn] = std::max(columnLower[iColumn], -dValue);
                columnUpper[iColumn] = std::min(columnUpper[iColumn], dValue);
             }
          }
      } break;
      case ClpParam::REALLY_SCALE:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          ClpSimplex newModel(model_, model_.scalingFlag());
          printGeneralMessage(model_, "model really really scaled\n");
          model_ = newModel;
        } break;
      case ClpParam::USERCLP:
        // Replace the sample code by whatever you want
        if (goodModel) {
          ClpSimplex *thisModel = &model_;
          buffer.str("");
          buffer << "Dummy user code - model has " << model_.numberRows()
                 << " rows and " << model_.numberColumns() << " columns"
                 << std::endl;
          printGeneralMessage(model_, buffer.str());
        }
        break;
      case ClpParam::HELP:
        std::cout << "Coin LP version " << CLP_VERSION << ", build " << __DATE__
                  << std::endl;
        std::cout << "Non default values:-" << std::endl;
        std::cout << "Perturbation " << model_.perturbation()
                  << " (default 100)" << std::endl;
        CoinParamUtils::printString("Presolve being done with 5 passes\n\
Dual steepest edge steep/partial on matrix shape and factorization density\n\
Clpnnnn taken out of messages\n\
If Factorization frequency default then done on size of matrix\n\n\
(-)unitTest, (-)netlib or (-)netlibp will do standard tests\n\n\
You can switch to interactive mode at any time so\n\
clp watson.mps -scaling off -primalsimplex\nis the same as\n\
clp watson.mps -\nscaling off\nprimalsimplex");
        break;
      case ClpParam::PRINTSOL:
      case ClpParam::WRITESOL:
      case ClpParam::WRITEGMPLSOL:{
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          fp = NULL;
          bool append = false;
          canOpen = false;
          if (paramCode == ClpParam::PRINTSOL){
             fp = stdout;
             fprintf(fp, "\n");
          } else {
             param->readValue(inputQueue, fileName, &message);
             if (fileName == "append"){
                switch (paramCode){
                 case ClpParam::WRITESOL:
                   fileName = parameters[ClpParam::SOLUTIONFILE]->fileName();
                   break;
                 case ClpParam::WRITEGMPLSOL:
                   fileName = parameters[ClpParam::GMPLSOLFILE]->fileName();
                   break;
                }
                CoinParamUtils::processFile(fileName,
                                  parameters[ClpParam::DIRECTORY]->dirName(),
                                            &canOpen);
                if (!canOpen) {
                   buffer.str("");
                   buffer << "Unable to open file " << fileName.c_str();
                   printGeneralMessage(model_, buffer.str());
                   continue;
                }
                append = true;
             } else if (paramCode == ClpParam::WRITESOL) {
                CoinParamUtils::processFile(fileName,
                                 parameters[ClpParam::DIRECTORY]->dirName());
                if (fileName == ""){
                   fileName = parameters[ClpParam::SOLUTIONFILE]->fileName();
                }else{
                   parameters[ClpParam::SOLUTIONFILE]->setFileName(fileName);
                }
             } else if (paramCode == ClpParam::WRITEGMPLSOL) {
                CoinParamUtils::processFile(fileName,
                                 parameters[ClpParam::DIRECTORY]->dirName());
                if (fileName == ""){
                   fileName = parameters[ClpParam::GMPLSOLFILE]->fileName();
                }else{
                   parameters[ClpParam::GMPLSOLFILE]->setFileName(fileName);
                }
             }

             if (!append) {
                fp = fopen(fileName.c_str(), "w");
             } else {
                fp = fopen(fileName.c_str(), "a");
             }
             if (!fp){
                buffer.str("");
                buffer << "Unable to open file " << fileName << std::endl;
                printGeneralWarning(model_, buffer.str());
             }
          }
            // See if Glpk
            if (paramCode == ClpParam::WRITEGMPLSOL) {
              int numberRows = model_.getNumRows();
              int numberColumns = model_.getNumCols();
              int numberGlpkRows = numberRows + 1;
              fprintf(fp, "%d %d\n", numberGlpkRows, numberColumns);
              int iStat = model_.status();
              int iStat2 = GLP_UNDEF;
              if (iStat == 0) {
                // optimal
                iStat2 = GLP_FEAS;
              } else if (iStat == 1) {
                // infeasible
                iStat2 = GLP_NOFEAS;
              } else if (iStat == 2) {
                // unbounded
                // leave as 1
              } else if (iStat >= 3 && iStat <= 5) {
                iStat2 = GLP_FEAS;
              }
              double objValue =
                  model_.getObjValue() * model_.getObjSense();
              fprintf(fp, "%d 2 %g\n", iStat2, objValue);
              if (numberGlpkRows > numberRows) {
                // objective as row
                fprintf(fp, "4 %g 1.0\n", objValue);
              }
              int lookup[6] = {4, 1, 3, 2, 4, 5};
              const double *primalRowSolution =
                  model_.primalRowSolution();
              const double *dualRowSolution = model_.dualRowSolution();
              for (int i = 0; i < numberRows; i++) {
                fprintf(fp, "%d %g %g\n",
                        lookup[model_.getRowStatus(i)],
                        primalRowSolution[i], dualRowSolution[i]);
              }
              const double *primalColumnSolution =
                  model_.primalColumnSolution();
              const double *dualColumnSolution =
                  model_.dualColumnSolution();
              for (int i = 0; i < numberColumns; i++) {
                fprintf(fp, "%d %g %g\n",
                        lookup[model_.getColumnStatus(i)],
                        primalColumnSolution[i], dualColumnSolution[i]);
              }
              fclose(fp);
              break;
            }
            // Write solution header (suggested by Luigi Poderico)
            double objValue = model_.getObjValue();
            int iStat = model_.status();
            if (iStat == 0) {
              fprintf(fp, "Optimal");
            } else if (iStat == 1) {
              // infeasible
              fprintf(fp, "Infeasible");
            } else if (iStat == 2) {
              // unbounded
              fprintf(fp, "Unbounded");
            } else if (iStat == 3) {
              fprintf(fp, "Stopped on iterations or time");
            } else if (iStat == 4) {
              fprintf(fp, "Stopped on difficulties");
            } else if (iStat == 5) {
              fprintf(fp, "Stopped on ctrl-c");
            } else {
              fprintf(fp, "Status unknown");
            }
            char printFormat[50];
            sprintf(printFormat, " - objective value %s\n",
                    CLP_QUOTE(CLP_OUTPUT_FORMAT));
            fprintf(fp, printFormat, objValue);
            if (printMode == 9) {
              // just statistics
              int numberRows = model_.numberRows();
              double *dualRowSolution = model_.dualRowSolution();
              double *primalRowSolution = model_.primalRowSolution();
              double *rowLower = model_.rowLower();
              double *rowUpper = model_.rowUpper();
              double highestPrimal;
              double lowestPrimal;
              double highestDual;
              double lowestDual;
              double largestAway;
              int numberAtLower;
              int numberAtUpper;
              int numberBetween;
              highestPrimal = -COIN_DBL_MAX;
              lowestPrimal = COIN_DBL_MAX;
              highestDual = -COIN_DBL_MAX;
              lowestDual = COIN_DBL_MAX;
              largestAway = 0.0;
              ;
              numberAtLower = 0;
              numberAtUpper = 0;
              numberBetween = 0;
              for (int iRow = 0; iRow < numberRows; iRow++) {
                double primal = primalRowSolution[iRow];
                double lower = rowLower[iRow];
                double upper = rowUpper[iRow];
                double dual = dualRowSolution[iRow];
                highestPrimal = std::max(highestPrimal, primal);
                lowestPrimal = std::min(lowestPrimal, primal);
                highestDual = std::max(highestDual, dual);
                lowestDual = std::min(lowestDual, dual);
                if (primal < lower + 1.0e-6) {
                  numberAtLower++;
                } else if (primal > upper - 1.0e-6) {
                  numberAtUpper++;
                } else {
                  numberBetween++;
                  largestAway = std::max(
                      largestAway, std::min(primal - lower, upper - primal));
                }
              }
              buffer.str("");
              buffer << "For rows " << numberAtLower << " at lower, "
                     << numberBetween << " between, " << numberAtUpper
                     << " at upper - lowest " << lowestPrimal << ", highest "
                     << highestPrimal << " most away " << largestAway
                     << " - highest dual " << highestDual << " lowest "
                     << lowestDual << std::endl;
              printGeneralMessage(model_, buffer.str());
              int numberColumns = model_.numberColumns();
              double *dualColumnSolution = model_.dualColumnSolution();
              double *primalColumnSolution =
                  model_.primalColumnSolution();
              double *columnLower = model_.columnLower();
              double *columnUpper = model_.columnUpper();
              highestPrimal = -COIN_DBL_MAX;
              lowestPrimal = COIN_DBL_MAX;
              highestDual = -COIN_DBL_MAX;
              lowestDual = COIN_DBL_MAX;
              largestAway = 0.0;
              ;
              numberAtLower = 0;
              numberAtUpper = 0;
              numberBetween = 0;
              for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                double primal = primalColumnSolution[iColumn];
                double lower = columnLower[iColumn];
                double upper = columnUpper[iColumn];
                double dual = dualColumnSolution[iColumn];
                highestPrimal = std::max(highestPrimal, primal);
                lowestPrimal = std::min(lowestPrimal, primal);
                highestDual = std::max(highestDual, dual);
                lowestDual = std::min(lowestDual, dual);
                if (primal < lower + 1.0e-6) {
                  numberAtLower++;
                } else if (primal > upper - 1.0e-6) {
                  numberAtUpper++;
                } else {
                  numberBetween++;
                  largestAway = std::max(
                      largestAway, std::min(primal - lower, upper - primal));
                }
              }
              buffer.str("");
              buffer << "For columns " << numberAtLower << " at lower, "
                     << numberBetween << " between, " << numberAtUpper
                     << " at upper - lowest " << lowestPrimal << ", highest "
                     << highestPrimal << " most away " << largestAway
                     << " - highest dual " << highestDual << " lowest "
                     << lowestDual << std::endl;
              printGeneralMessage(model_, buffer.str());
              break;
            }
            // make fancy later on
            int iRow;
            int numberRows = model_.numberRows();
            int lengthName = model_.lengthNames(); // 0 if no names
            int lengthPrint = std::max(lengthName, 8);
            // in general I don't want to pass around massive
            // amounts of data but seems simpler here
            std::vector<std::string> rowNames = *(model_.rowNames());
            std::vector<std::string> columnNames =
                *(model_.columnNames());

            double *dualRowSolution = model_.dualRowSolution();
            double *primalRowSolution = model_.primalRowSolution();
            double *rowLower = model_.rowLower();
            double *rowUpper = model_.rowUpper();
            double primalTolerance = model_.primalTolerance();
            bool doMask = (parameters[ClpParam::PRINTMASK]->strVal() != "" &&
                           lengthName);
            int *maskStarts = NULL;
            int maxMasks = 0;
            char **masks = NULL;
            if (doMask) {
              int nAst = 0;
              const char *pMask2 =
                 parameters[ClpParam::PRINTMASK]->strVal().c_str();
              char pMask[100];
              size_t lengthMask = strlen(pMask2);
              assert(lengthMask < 100);
              if (*pMask2 == '"') {
                if (pMask2[lengthMask - 1] != '"') {
                  buffer.str("");
                  buffer << "mismatched \" in mask " << pMask2 << std::endl;
                  printGeneralWarning(model_, buffer.str());
                  break;
                } else {
                  strcpy(pMask, pMask2 + 1);
                  *strchr(pMask, '"') = '\0';
                }
              } else if (*pMask2 == '\'') {
                if (pMask2[lengthMask - 1] != '\'') {
                  buffer.str("");
                  buffer << "mismatched ' in mask " << pMask2 << std::endl;
                  printGeneralWarning(model_, buffer.str());
                  break;
                } else {
                  strcpy(pMask, pMask2 + 1);
                  *strchr(pMask, '\'') = '\0';
                }
              } else {
                strcpy(pMask, pMask2);
              }
              if (lengthMask > static_cast<size_t>(lengthName)) {
                buffer.str("");
                buffer << "mask " << pMask << " too long - skipping"
                       << std::endl;
                printGeneralWarning(model_, buffer.str());
                break;
              }
              maxMasks = 1;
              for (size_t iChar = 0; iChar < lengthMask; iChar++) {
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
                  for (size_t i = 0; i <= lengthName - length; i++) {
                    char *maskOut = newMasks[nEntries];
                    CoinMemcpyN(oldMask, static_cast<int>(nBefore), maskOut);
                    for (size_t k = 0; k < i; k++)
                      maskOut[k + nBefore] = '?';
                    CoinMemcpyN(ast + 1, static_cast<int>(nAfter),
                                maskOut + nBefore + i);
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
            if (printMode > 5) {
              int numberColumns = model_.numberColumns();
              // column length unless rhs ranging
              int number = numberColumns;
              switch (printMode) {
              // bound ranging
              case 6:
                fprintf(fp, "Bound ranging");
                break;
              // rhs ranging
              case 7:
                fprintf(fp, "Rhs ranging");
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
                    printGeneralWarning(model_, "No names match - doing all\n");
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
                    printGeneralWarning(model_, "No names match - doing all\n");
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
                model_.primalRanging(numberRows, which, valueIncrease,
                                             sequenceIncrease, valueDecrease,
                                             sequenceDecrease);
                break;
              // objective ranging
              case 8:
                model_.dualRanging(number, which, valueIncrease,
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
            sprintf(printFormat, " %s         %s\n",
                    CLP_QUOTE(CLP_OUTPUT_FORMAT), CLP_QUOTE(CLP_OUTPUT_FORMAT));
            if (printMode > 2) {
              for (iRow = 0; iRow < numberRows; iRow++) {
                int type = printMode - 3;
                if (primalRowSolution[iRow] >
                        rowUpper[iRow] + primalTolerance ||
                    primalRowSolution[iRow] <
                        rowLower[iRow] - primalTolerance) {
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
                    const char *name = rowNames[iRow].c_str();
                    size_t n = strlen(name);
                    size_t i;
                    for (i = 0; i < n; i++)
                      fprintf(fp, "%c", name[i]);
                    for (; i < static_cast<size_t>(lengthPrint); i++)
                      fprintf(fp, " ");
                  }
                  fprintf(fp, printFormat, primalRowSolution[iRow],
                          dualRowSolution[iRow]);
                }
              }
            }
            int iColumn;
            int numberColumns = model_.numberColumns();
            double *dualColumnSolution = model_.dualColumnSolution();
            double *primalColumnSolution =
                model_.primalColumnSolution();
            double *columnLower = model_.columnLower();
            double *columnUpper = model_.columnUpper();
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
              int type = (printMode > 3) ? 1 : 0;
              if (primalColumnSolution[iColumn] >
                      columnUpper[iColumn] + primalTolerance ||
                  primalColumnSolution[iColumn] <
                      columnLower[iColumn] - primalTolerance) {
                fprintf(fp, "** ");
                type = 2;
              } else if (fabs(primalColumnSolution[iColumn]) > 1.0e-8) {
                type = 1;
              } else if (numberColumns < 50) {
                type = 3;
              }
              if (doMask &&
                  !maskMatches(maskStarts, masks, columnNames[iColumn]))
                type = 0;
              if (type) {
                fprintf(fp, "%7d ", iColumn);
                if (lengthName) {
                  const char *name = columnNames[iColumn].c_str();
                  size_t n = strlen(name);
                  size_t i;
                  for (i = 0; i < n; i++)
                    fprintf(fp, "%c", name[i]);
                  for (; i < static_cast<size_t>(lengthPrint); i++)
                    fprintf(fp, " ");
                }
                fprintf(fp, printFormat, primalColumnSolution[iColumn],
                        dualColumnSolution[iColumn]);
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
      case ClpParam::WRITESOLBINARY: {
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          // get next field
           param->readValue(inputQueue, fileName, &message);
           CoinParamUtils::processFile(fileName,
                                 parameters[ClpParam::DIRECTORY]->dirName());
           if (fileName == ""){
              fileName = parameters[ClpParam::SOLUTIONBINARYFILE]->fileName();
           }else{
              parameters[ClpParam::SOLUTIONBINARYFILE]->setFileName(fileName);
           }
           ClpParamUtils::saveSolution(&model_, fileName);
        } break;
      case ClpParam::ENVIRONMENT: {
#if !defined(_MSC_VER) && !defined(__MSVCRT__)
        // Don't think it will work with Visual Studio
        char *input = getenv("CLP_ENVIRONMENT");
        if (input) {
          while (!inputQueue.empty()){
             inputQueue.pop_front();
          }
          std::istringstream inputStream(input);
          CoinParamUtils::readFromStream(inputQueue, inputStream);
        }
#else
        printGeneralWarning(
            model_, "** Parameter not valid on Windows with Visual Studio\n");
#endif
        break;
      }
      case ClpParam::GUESS: {
        if (!goodModel){
          printGeneralWarning(model_, "** Current model not valid\n");
          continue;
        }
          ClpSimplexOther *model2 =
              static_cast<ClpSimplexOther *>(&model_);
          std::string input = model2->guess(0);
          if (input != "") {
            while (!inputQueue.empty()){
               inputQueue.pop_front();
            }
            std::istringstream inputStream(input);
            CoinParamUtils::readFromStream(inputQueue, inputStream);
          } else {
            printGeneralWarning(model_,
                                "** Guess unable to generate commands\n");
          }

      }  break;
      default:
        abort();
      }
    }
  }
  // By now all memory should be freed
#ifdef DMALLOC
  dmalloc_log_unfreed();
  dmalloc_shutdown();
#endif
#ifdef CLP_USEFUL_PRINTOUT
  printf(
      "BENCHMARK_RESULT %s took %g cpu seconds (%g elapsed - %d iterations) "
      "- %d row, %d columns, %d elements - reduced to %d, %d, %d\n",
      mpsFile.c_str(), CoinCpuTime() - startCpu,
      CoinGetTimeOfDay() - startElapsed, debugInt[23], debugInt[0],
      debugInt[1], debugInt[2], debugInt[3], debugInt[4], debugInt[5]);
  const char *method[] = {"", "Dual", "Primal", "Sprint", "Idiot"};
  printf("BENCHMARK using method %s (%d)\n", method[debugInt[6]],
         debugInt[7]);
#endif
  return 0;
}

int clpReadAmpl(ampl_info * info, int argc, char **argv, ClpSimplex &model)
{

   std::ostringstream buffer;
   memset(&info, 0, sizeof(info));
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

   union {
      void *voidModel;
      CoinModel *model;
   } coinModelStart;
   coinModelStart.model = NULL;
   int returnCode =
      readAmpl(info, argc, argv, &coinModelStart.voidModel, "clp");
   if (returnCode) {
      return returnCode;
   }
   if (info->numberBinary + info->numberIntegers + info->numberSos &&
       !info->starts) {
      printGeneralWarning(model, "Unable to handle integer problems\n");
      return 1;
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
    if (!noPrinting) {
      buffer.str("");
      buffer << info->numberRows << " rows, " << info->numberColumns
             << " columns and " << info->numberElements << " elements"
             << std::endl;
      printGeneralMessage(model, buffer.str());
    }
    if (!coinModelStart.model) {
      // linear
      model.loadProblem(info->numberColumns, info->numberRows,
                          reinterpret_cast<const CoinBigIndex *>(info->starts),
                          info->rows, info->elements, info->columnLower,
                          info->columnUpper, info->objective, info->rowLower,
                          info->rowUpper);
    } else {
      // QP
      model.loadProblem(*(coinModelStart.model));
    }
    // If we had a solution use it
    if (info->primalSolution) {
      model.setColSolution(info->primalSolution);
    }
    // status
    if (info->rowStatus) {
      unsigned char *statusArray = model.statusArray();
      int i;
      for (i = 0; i < info->numberColumns; i++)
        statusArray[i] = static_cast<unsigned char>(info->columnStatus[i]);
      statusArray += info->numberColumns;
      for (i = 0; i < info->numberRows; i++)
        statusArray[i] = static_cast<unsigned char>(info->rowStatus[i]);
    }
    freeArrays1(info);
    // modify objective if necessary
    model.setOptimizationDirection(info->direction);
    model.setObjectiveOffset(-info->offset);
    if (info->offset) {
      buffer.str("");
      buffer << "Ampl objective offset is " << info->offset << std::endl;
      printGeneralMessage(model, buffer.str());
    }
    return 0;
  }

  static void breakdown(const char *name, CoinBigIndex numberLook,
                        const double *region) {
    double range[] = {-COIN_DBL_MAX, -1.0e15,     -1.0e11,  -1.0e8,  -1.0e5,
                      -1.0e4,        -1.0e3,      -1.0e2,   -1.0e1,  -1.0,
                      -1.0e-1,       -1.0e-2,     -1.0e-3,  -1.0e-4, -1.0e-5,
                      -1.0e-8,       -1.0e-11,    -1.0e-15, 0.0,     1.0e-15,
                      1.0e-11,       1.0e-8,      1.0e-5,   1.0e-4,  1.0e-3,
                      1.0e-2,        1.0e-1,      1.0,      1.0e1,   1.0e2,
                      1.0e3,         1.0e4,       1.0e5,    1.0e8,   1.0e11,
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

void sortOnOther(int *column, const CoinBigIndex *rowStart, int *order,
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

static void statistics(ClpSimplex * originalModel, ClpSimplex * model) {
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
              int *columnAdd = new int[spaceNeeded];
              CoinBigIndex *startAdd = new CoinBigIndex[numberDifferentObj + 1];
              double *elementAdd = new double[spaceNeeded];
              CoinBigIndex *rowAdd =
                  new CoinBigIndex[2 * numberDifferentObj + 1];
              int *newIsInteger =
                  reinterpret_cast<int *>(rowAdd + numberDifferentObj + 1);
              double *objectiveNew = new double[3 * numberDifferentObj];
              double *lowerNew = objectiveNew + numberDifferentObj;
              double *upperNew = lowerNew + numberDifferentObj;
              memset(startAdd, 0,
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
                      if (objectiveNew[numberDifferentObj] *
                              objective[iColumn] <
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
                                   objectiveNew, startAdd, NULL, NULL);
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
              delete[] startAdd;
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
          0 0/inf, 1 0/up, 2 lo/inf, 3 lo/up, 4 free, 5 fix, 6 -inf/0, 7
       -inf/up, 8 0/1
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
    int kMax = model->logLevel() > 3 ? 10000000 : 10;
    k = 0;
    for (iRow = 1; iRow <= numberRows; iRow++) {
      if (number[iRow]) {
        k++;
        printf("%d columns have %d entries\n", number[iRow], iRow);
        if (k == kMax)
          break;
      }
    }
    {
      int n = 0;
      int nLast = -1;
      int iLast = -1;
      int jRow = iRow;
      iRow++;
      for (; iRow < numberRows; iRow++) {
        if (number[iRow]) {
          n += number[iRow];
          nLast = number[iRow];
          iLast = iRow;
        }
      }
      if (n) {
        printf("\n    ... set logLevel >3 to see all ......\n\n");
        printf("%d columns > %d entries < %d\n", n - nLast, jRow, iLast);
        printf("%d column%s %d entries\n", nLast, nLast > 1 ? "s have" : " has",
               iLast);
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
    int logLevel = model->logLevel();
    bool reduceMaster = (logLevel & 16) != 0;
    bool writeMatrices = (logLevel & 8) != 0;
    bool morePrint = (logLevel & 32) != 0;
    if (logLevel > 3) {
      // See what is needed
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
            printf("Singleton free %d - cost %g\n", iColumn,
                   objective[iColumn]);
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
                   (rowLower[iRow0] == rowUpper[iRow0]) ? "==" : "", n0,
                   element1, (rowLower[iRow1] == rowUpper[iRow1]) ? "==" : "",
                   n1);
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
            if (morePrint)
              printf("Singleton %d with no objective in row with %d elements - "
                     "rhs %g,%g\n",
                     iColumn, rowLength[iRow], rowLower[iRow], rowUpper[iRow]);
            nPossibleZeroCost++;
          } else if (value != -COIN_DBL_MAX) {
            if (morePrint)
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
      int *blockStart = new int[3 * (numberRows + numberColumns) + 1];
      int *columnBlock = blockStart + numberRows;
      int *nextColumn = columnBlock + numberColumns;
      int *blockCount = nextColumn + numberColumns;
      int *blockEls = blockCount + numberRows + 1;
      int *countIntegers = blockEls + numberRows;
      memset(countIntegers, 0, numberColumns * sizeof(int));
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
              assert(columnBlock[iColumn] < 0 ||
                     columnBlock[iColumn] == iBlock);
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
        // If wanted minimize master rows
        if (reduceMaster)
          thisBestValue =
              (increment < 0) ? thisBestBreak : numberRows - thisBestBreak;
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
            firstMaster = thisBestBreak + 1;
            lastMaster = numberRows;
          }
        }
      }
      if (firstMaster <= lastMaster) {
        if (firstMaster == lastMaster)
          printf("Total decomposition! - ");
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
        int numberMasterIntegers = 0;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          int iBlock = columnBlock[iColumn];
          bool isInteger = (model->isInteger(iColumn));
          if (iBlock >= 0) {
            nextColumn[iBlock]++;
            if (isInteger)
              countIntegers[iBlock]++;
          } else {
            if (isInteger)
              numberMasterIntegers++;
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
               "%d, %d master columns (%d empty, %d integer) out of %d\n",
               useful ? "**Useful" : "NoGood", numberBlocks, largestRows,
               largestColumns, numberMaster, numberEmpty, numberRows,
               numberMasterColumns, numberEmptyColumns, numberMasterIntegers,
               numberColumns);
        for (int i = 0; i < numberBlocks; i++)
          printf("Block %d has %d rows and %d columns (%d elements, %d "
                 "integers)\n",
                 i, blockCount[i], nextColumn[i], blockEls[i],
                 countIntegers[i]);
        FILE *fpBlocks = fopen("blocks.data", "wb");
        printf("Blocks data on file blocks.data\n");
        int stats[3];
        stats[0] = numberRows;
        stats[1] = numberColumns;
        stats[2] = numberBlocks;
        size_t numberWritten;
        numberWritten = fwrite(stats, sizeof(int), 3, fpBlocks);
        assert(numberWritten == 3);
        numberWritten = fwrite(blockStart, sizeof(int), numberRows, fpBlocks);
        assert(numberWritten == numberRows);
        numberWritten =
            fwrite(columnBlock, sizeof(int), numberColumns, fpBlocks);
        assert(numberWritten == numberColumns);
        fclose(fpBlocks);
        if (writeMatrices) {
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
            for (int jColumn = 0; jColumn < nColumns; jColumn++) {
              int iColumn = whichColumns[jColumn];
              if (model->isInteger(iColumn))
                subset.setInteger(jColumn);
              std::string name = model->getColumnName(iColumn);
              subset.setColumnName(jColumn, name);
            }
            subset.writeMps(name, 0, 1);
          }
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
    {
      int n = 0;
      int nLast = -1;
      int iLast = -1;
      int jColumn = iColumn;
      iColumn++;
      for (; iColumn < numberColumns; iColumn++) {
        if (number[iColumn]) {
          n += number[iColumn];
          nLast = number[iColumn];
          iLast = iColumn;
        }
      }
      if (n) {
        printf("\n    ... set logLevel >3 to see all ......\n\n");
        printf("%d rows > %d entries < %d\n", n - nLast, jColumn, iLast);
        printf("%d row%s %d entries\n", nLast, nLast > 1 ? "s have" : " has",
               iLast);
      }
    }
    if (morePrint
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
              if (morePrint) {
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
            for (int jLook = iLook + 1; jLook < iLook + numberPossible;
                 jLook++) {
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
            printf("possible %d - %d integers\n", numberPossible,
                   numberIntegers);
            int lastLook = iLook;
            int nMapRow = -1;
            for (int jLook = iLook + 1; jLook < iLook + numberPossible;
                 jLook++) {
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
                    assert(iRow == mapRow[row2]);
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
                      if (model->isInteger(iColumn) !=
                          model->isInteger(iColumn2))
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
        printf("Adding %d rows with %d elements\n", nAddRows,
               startAdd[nAddRows]);
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
    breakdown("Elements", numberElements, elementByColumn);
    breakdown("RowLower", numberRows, rowLower);
    breakdown("RowUpper", numberRows, rowUpper);
    breakdown("ColumnLower", numberColumns, columnLower);
    breakdown("ColumnUpper", numberColumns, columnUpper);
    breakdown("Objective", numberColumns, objective);
    // do integer objective
    double *obj = CoinCopyOfArray(objective, numberColumns);
    int n = 0;
    //#define FIX_COSTS 1.0
#ifdef FIX_COSTS
    double *obj2 = originalModel->objective();
    double *upper2 = originalModel->columnUpper();
#endif
    for (int i = 0; i < numberColumns; i++) {
      if (integerInformation && integerInformation[i]) {
        obj[n++] = obj[i];
#ifdef FIX_COSTS
        if (obj2[i] > FIX_COSTS)
          upper2[i] = 0.0;
#endif
      }
    }
    if (n) {
      breakdown("Integer objective", n, obj);
    }
}

static bool maskMatches(const int *starts, char **masks, std::string &check) {
    // back to char as I am old fashioned
    const char *checkC = check.c_str();
    size_t length = strlen(checkC);
    while (checkC[length - 1] == ' ')
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

static void clean(char *temp) {
    char *put = temp;
    while (*put >= ' ')
      put++;
    *put = '\0';
}

static void generateCode(const char *fileName, int type) {
    FILE *fp = fopen(fileName, "r");
    assert(fp);
    int numberLines = 0;
#define MAXLINES 500
#define MAXONELINE 200
    char line[MAXLINES][MAXONELINE];
    while (fgets(line[numberLines], MAXONELINE, fp)) {
      assert(numberLines < MAXLINES);
      clean(line[numberLines]);
      numberLines++;
    }
    fclose(fp);
    // add in actual solve
    strcpy(line[numberLines], "5  clpModel->initialSolve(clpSolve);");
    numberLines++;
    fp = fopen(fileName, "w");
    assert(fp);
    char apo = '"';
    char backslash = '\\';

    fprintf(fp, "#include %cClpSimplex.hpp%c\n", apo, apo);
    fprintf(fp, "#include %cClpSolve.hpp%c\n", apo, apo);

    fprintf(fp, "\nint main (int argc, const char *argv[])\n{\n");
    fprintf(fp, "  ClpSimplex  model;\n");
    fprintf(fp, "  int status=1;\n");
    fprintf(fp, "  if (argc<2)\n");
    fprintf(fp, "    fprintf(stderr,%cPlease give file name%cn%c);\n", apo,
            backslash, apo);
    fprintf(fp, "  else\n");
    fprintf(fp, "    status=model.readMps(argv[1],true);\n");
    fprintf(fp, "  if (status) {\n");
    fprintf(fp, "    fprintf(stderr,%cBad readMps %%s%cn%c,argv[1]);\n", apo,
            backslash, apo);
    fprintf(fp, "    exit(1);\n");
    fprintf(fp, "  }\n\n");
    fprintf(fp, "  // Now do requested saves and modifications\n");
    fprintf(fp, "  ClpSimplex * clpModel = & model;\n");
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
                             "Add to model"};
    for (int iType = 0; iType < 9; iType++) {
      if (!wanted[iType])
        continue;
      int n = 0;
      int iLine;
      for (iLine = 0; iLine < numberLines; iLine++) {
        if (line[iLine][0] == '0' + iType) {
          if (!n)
            fprintf(fp, "\n  // %s\n\n", header[iType].c_str());
          n++;
          fprintf(fp, "%s\n", line[iLine] + 1);
        }
      }
    }
    fprintf(fp, "\n  // Now you would use solution etc etc\n\n");
    fprintf(fp, "  return 0;\n}\n");
    fclose(fp);
    printf("C++ file written to %s\n", fileName);
}

  /* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
   */
