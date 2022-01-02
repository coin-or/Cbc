/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <cassert>
#include <string>
#include <sstream>

#include "CoinUtilsConfig.h"

#include "CoinParam.hpp"
#include "CoinFileIO.hpp"
#include "CoinFinite.hpp"

#include "CbcModel.hpp"

#include "CbcParam.hpp"
#include "CbcParamUtils.hpp"
#include "CbcParameters.hpp"

namespace CbcParamUtils {

//###########################################################################
//###########################################################################

/* Functions to perform actions related to setting parameters */

/*
  Maintainer's utility, scan the parameters and report the ones that are
  unimplemented (i.e., have no pushFunc).
*/

int doUnimplementedParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters parameters = *cbcParam.parameters();

  int unimpCnt = 0;
  int maxAcross = 5;
  for (int code = CbcParam::FIRSTPARAM + 1;
       code < CbcParam::LASTPARAM; code++) {
    if (parameters[code]->pushFunc() == 0) {
      if (unimpCnt % maxAcross == 0) {
        std::cout << std::endl;
      } else {
        std::cout << " ";
      }
      std::cout << parameters[code]->name();
      unimpCnt++;
    }
  }
  if (unimpCnt % maxAcross != 1) {
    std::cout << std::endl;
  }
  std::cout << unimpCnt << " unimplemented parameters." << std::endl;

  return (0);
}

//###########################################################################
//###########################################################################

/*
  Noop function. Mainly to eliminate commands from the list returned by
  doUnimplmentedParam.
*/

int doNothingParam(CoinParam &param) { return (0); }

/*
  Function to terminate command parsing by returning -1.
*/

int doExitParam(CoinParam &param)

{
  return (-1);
}

//###########################################################################
//###########################################################################

/*
  Function to print the current version.
*/

int doVersionParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  std::cout << "Cbc version " << parameters->getVersion() << std::endl;

  return (0);
}

//###########################################################################
//###########################################################################

/*
  Function to handle help (HELP), `?' (GENERALQUERY), and `???'
  (FULLGENERALQUERY).
*/

int doHelpParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  int cbcParamCode = cbcParam.paramCode();

  int verbose = parameters->getVerbose();
  bool shortHelp = ((verbose & 0x01) ? true : false);
  bool longHelp = ((verbose & 0x02) ? true : false);
  bool hidden = ((verbose & 0x08) ? true : false);

  CoinParamVec &paramVec = parameters->paramVec();
  /*
      Tune up the initial parameters. FULLGENERALQUERY will print normally
     hidden params, and a request for long help overrules a request for short
     help.
    */
  if (cbcParamCode == CbcParam::FULLGENERALQUERY) {
    hidden = true;
  }
  if (longHelp) {
    shortHelp = false;
  }

  CoinParamUtils::printGenericHelp();

  std::cout << "\nAvailable commands are:";
  std::string pfx("  ");
  CoinParamUtils::printHelp(paramVec, 0, paramVec.size() - 1,
                            pfx, shortHelp, longHelp, hidden);

  return (0);
}

//###########################################################################
//###########################################################################

/*
  Function to push a double-valued parameter.
*/

int pushCbcSolverDblParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  double val = cbcParam.dblVal();
  int cbcParamCode = cbcParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (cbcParamCode) {
   case CbcParam::ARTIFICIALCOST: {
      parameters->setArtVarMode(CbcParameters::ParamOn, val);
      break;
   }
   case CbcParam::DEXTRA3: {
      parameters->setExtraDbl3(val);
      break;
   }
   case CbcParam::DEXTRA4: {
      parameters->setExtraDbl4(val);
      break;
   }
   case CbcParam::DEXTRA5: {
      parameters->setExtraDbl5(val);
      break;
   }
   case CbcParam::DJFIX: {
      parameters->setDjFixMode(CbcParameters::ParamOn, val);
      break;
   }
   case CbcParam::FAKECUTOFF: {
      parameters->setFeasPumpFakeCutoff(val);
      break;
   }
   case CbcParam::FAKEINCREMENT: {
      parameters->setFeasPumpFakeIncrement(val);
      break;
   }
   case CbcParam::SMALLBAB: {
      parameters->setSmallBab(val);
      break;
   }
   case CbcParam::TIGHTENFACTOR: {
      parameters->setTightenFactor(val);
      break;
   }
     
   default: {
      std::cerr << "pushCbcSolverDbl: no equivalent CbcParameters field for "
                << "parameter code `" << cbcParamCode << "'." << std::endl;
      retval = -1;
      break;
   }
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push an integer-valued parameter.
*/

int pushCbcSolverIntParam(CoinParam &param)
{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  int val = cbcParam.intVal();
  int cbcParamCode = cbcParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (cbcParamCode) {
  case CbcParam::BKPIVOTINGSTRATEGY: {
    parameters->setBkPivotStrategy(val);
    break;
  }
  case CbcParam::BKMAXCALLS: {
    parameters->setBkMaxCalls(val);
    break;
  }
  case CbcParam::BKCLQEXTMETHOD: {
    parameters->setBkClqExtMethod(val);
    break;
  }
  case CbcParam::CPP: {
    parameters->setCppMode(val);
    break;
  }
  case CbcParam::CUTDEPTH: {
    parameters->setCutDepth(val);
    break;
  }
  case CbcParam::CUTLENGTH: {
    parameters->setCutLength(val);
    break;
  }
  case CbcParam::CUTPASSINTREE: {
    parameters->setCutPassInTree(val);
    break;
  }
  case CbcParam::DEPTHMINIBAB: {
    parameters->setDepthMiniBaB(val);
    break;
  }
  case CbcParam::DIVEOPT: {
    parameters->setDiveOpt(val);
    break;
  }
  case CbcParam::DIVEOPTSOLVES: {
    parameters->setDiveOptSolves(val);
    break;
  }
  case CbcParam::EXPERIMENT: {
    parameters->setExperimentMode(val);
    break;
  }
  case CbcParam::EXTRA1: {
    parameters->setExtraIntParam1(val);
    break;
  }
  case CbcParam::EXTRA2: {
    parameters->setExtraIntParam2(val);
    break;
  }
  case CbcParam::EXTRA3: {
    parameters->setExtraIntParam3(val);
    break;
  }
  case CbcParam::EXTRA4: {
    parameters->setExtraIntParam4(val);
    break;
  }
  case CbcParam::FPUMPITS: {
    parameters->setFeasPumpIters(val);
    break;
  }
  case CbcParam::FPUMPTUNE: {
    parameters->setFeasPumpTune(val);
    break;
  }
  case CbcParam::FPUMPTUNE2: {
    parameters->setFeasPumpTune2(val);
    break;
  }
  case CbcParam::HEUROPTIONS: {
    parameters->setHeurOptions(val);
    break;
  }
  case CbcParam::LOGLEVEL: {
    parameters->setLogLevel(val);
    break;
  }
  case CbcParam::LPLOGLEVEL: {
    parameters->setLpLogLevel(val);
    break;
  }
  case CbcParam::MAXSAVEDSOLS: {
    parameters->setMaxSavedSols(val);
    break;
  }
  case CbcParam::MAXSLOWCUTS: {
    parameters->setMaxSlowCuts(val);
    break;
  }
  case CbcParam::MOREMOREMIPOPTIONS: {
    parameters->setMoreMoreOptions(val);
    break;
  }
  case CbcParam::MULTIPLEROOTS: {
    parameters->setMultipleRoots(val);
    break;
  }
  case CbcParam::ODDWEXTMETHOD: {
    parameters->setOddWextMethod(val);
    break;
  }
  case CbcParam::OUTPUTFORMAT: {
    parameters->setOutputFormat(val);
    break;
  }
  case CbcParam::PRINTOPTIONS: {
    parameters->setPrintOptions(val);
    break;
  }
  case CbcParam::PROCESSTUNE: {
    parameters->setProcessTune(val);
    break;
  }
  case CbcParam::RANDOMSEED: {
    parameters->setRandomSeed(val);
    break;
  }
  case CbcParam::STRONGSTRATEGY: {
    parameters->setStrongStrategy(val);
    break;
  }
  case CbcParam::TESTOSI: {
    parameters->setTestOsi(val);
    break;
  }
#ifdef CBC_THREADS
  case CbcParam::THREADS: {
    parameters->setThreads(val);
    break;
  }
#endif
  case CbcParam::USERCBC: {
    parameters->setUserCbc(val);
    break;
  }
  case CbcParam::VERBOSE: {
    parameters->setVerbose(val);
    break;
  }
  case CbcParam::VUBTRY: {
    parameters->setVubTry(val);
    break;
  }
  default: {
    std::cerr << "pushCbcSolverInt: no method for storing "
              << "parameter code `" << cbcParamCode << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push a keyword-valued parameter. This is the catch-all function
  for keyword parameters that don't belong to any other useful grouping.
*/

int pushCbcSolverKwdParam(CoinParam &param) {
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  int mode = cbcParam.modeVal();
  int cbcParamCode = cbcParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
  */

  switch (cbcParamCode) {
  case CbcParam::COMMANDPRINTLEVEL: {
    parameters->setCommandDisplayMode(
                 static_cast<CbcParameters::CommandDisplayMode>(mode));
    break;
  }
  case CbcParam::CLQSTRENGTHENING: {
    parameters->setClqStrMode(static_cast<CbcParameters::ClqStrMode>(mode));
    break;
  }
  case CbcParam::BRANCHPRIORITY: {
    parameters->setBranchingPriorityMode(CbcParameters::BPOff);
    break;
  }
  case CbcParam::CUTOFFCONSTRAINT: {
    parameters->setCutoffMode(static_cast<CbcParameters::CutoffMode>(mode));
    break;
  }
  case CbcParam::INTPRINT: {
    parameters->setIntPrintMode(
        static_cast<CbcParameters::IntPrintMode>(mode));
    break;
  }
  case CbcParam::NODESTRATEGY: {
    parameters->setNodeStrategy(
        static_cast<CbcParameters::NodeStrategy>(mode));
    break;
  }
  case CbcParam::ORBITAL: {
    parameters->setOrbitalStrategy(
        static_cast<CbcParameters::OrbitalStrategy>(mode));
    break;
  }
  case CbcParam::PREPROCESS: {
    parameters->setIPPMode(static_cast<CbcParameters::IPPMode>(mode));
    break;
  }
  case CbcParam::SOSPRIORITIZE: {
    parameters->setSOSStrategy(static_cast<CbcParameters::SOSStrategy>(mode));
    break;
  }
  case CbcParam::STRATEGY: {
    parameters->setStrategyMode(
        static_cast<CbcParameters::StrategyMode>(mode));
    break;
  }
  case CbcParam::TIMEMODE: {
    parameters->setClockType(static_cast<CbcParameters::ClockType>(mode));
    break;
  }
  case CbcParam::USECGRAPH: {
    parameters->setCGraphMode(static_cast<CbcParameters::CGraphMode>(mode));
    break;
  }
  default:
    break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push a bool-valued parameter. These are really just keyword
  parameters that take values "on" and "off"

*/

int pushCbcSolverBoolParam(CoinParam &param) {
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  // This is ugly, get keyword and set parameter with it instead.
  CbcParameters::OnOffMode mode =
      static_cast<CbcParameters::OnOffMode>(cbcParam.modeVal());
  int cbcParamCode = cbcParam.paramCode();

  int retval = 0;

  switch (cbcParamCode) {
  case CbcParam::CPX: {
    parameters->setCPXMode(mode);
    break;
  }
  case CbcParam::DOHEURISTIC: {
    parameters->setDoHeuristicMode(mode);
    break;
  }
  case CbcParam::ERRORSALLOWED: {
    parameters->setImportErrorsMode(mode);
    break;
  }
  case CbcParam::MESSAGES: {
    parameters->setMessagePrefixMode(mode);
    break;
  }
  case CbcParam::PREPROCNAMES: {
    parameters->setPreProcNamesMode(mode);
    break;
  }
  case CbcParam::SOS: {
    parameters->setSOSMode(mode);
    break;
  }
  case CbcParam::USESOLUTION: {
    parameters->setUseSolutionMode(mode);
    break;
  }
  default:
    break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push a keyword-valued parameter related to heuristics.
*/

int pushCbcSolverHeurParam(CoinParam &param) {
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  // This is ugly, get keyword and set parameter with it instead.
  CbcParameters::HeurMode mode =
      static_cast<CbcParameters::HeurMode>(cbcParam.modeVal());
  int cbcParamCode = cbcParam.paramCode();

  int retval = 0;

  /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
  */
  switch (cbcParamCode) {
  case CbcParam::HEURISTICSTRATEGY: {
    parameters->setGreedyCoverMode(mode);
    parameters->setGreedyEqualityMode(mode);
    parameters->setCrossoverMode(mode);
    parameters->setDinsMode(mode);
    parameters->setDiveCoefficientMode(mode);
    parameters->setDiveFractionalMode(mode);
    parameters->setDiveGuidedMode(mode);
    parameters->setDiveLineSearchMode(mode);
    parameters->setDivePseudocostMode(mode);
    parameters->setDiveVectorLengthMode(mode);
    parameters->setDWMode(mode);
    parameters->setNaiveHeurMode(mode);
    parameters->setPivotAndFixMode(mode);
#if 0
      parameters->setPivotAndComplementMode(mode);
#endif
    parameters->setProximityMode(mode);
    parameters->setRandRoundMode(mode);
    parameters->setRensMode(mode);
    parameters->setRinsMode(mode);
    parameters->setRoundingMode(mode);
    parameters->setVndMode(mode);
    break;
  }
  case CbcParam::COMBINE: {
    parameters->setCombineMode(mode);
    break;
  }
  case CbcParam::CROSSOVER: {
    parameters->setCrossoverMode(mode);
    break;
  }
  case CbcParam::DINS: {
    parameters->setDinsMode(mode);
    break;
  }
  case CbcParam::DIVINGC: {
    parameters->setDiveCoefficientMode(mode);
    break;
  }
  case CbcParam::DIVINGF: {
    parameters->setDiveFractionalMode(mode);
    break;
  }
  case CbcParam::DIVINGG: {
    parameters->setDiveGuidedMode(mode);
    break;
  }
  case CbcParam::DIVINGL: {
    parameters->setDiveLineSearchMode(mode);
    break;
  }
  case CbcParam::DIVINGP: {
    parameters->setDivePseudocostMode(mode);
    break;
  }
  case CbcParam::DIVINGS: {
    parameters->setDiveRandomMode(mode);
    break;
  }
  case CbcParam::DIVINGV: {
    parameters->setDiveVectorLengthMode(mode);
    break;
  }
  case CbcParam::DW: {
    parameters->setDWMode(mode);
    break;
  }
  case CbcParam::FPUMP: {
    parameters->setFeasPumpMode(mode);
    break;
  }
  case CbcParam::GREEDY: {
    parameters->setGreedyCoverMode(mode);
    parameters->setGreedyEqualityMode(mode);
    break;
  }
  case CbcParam::NAIVE: {
    parameters->setNaiveHeurMode(mode);
    break;
  }
  case CbcParam::PIVOTANDFIX: {
    parameters->setPivotAndFixMode(mode);
    break;
  }
#if 0
   case CbcParam::PIVOTANDCOMPLEMENT: {
      parameters->setPivotAndComplementMode(mode);
      break;
   }
#endif
  case CbcParam::PROXIMITY: {
    parameters->setProximityMode(mode);
    break;
  }
  case CbcParam::RANDROUND: {
    parameters->setRandRoundMode(mode);
    break;
  }
  case CbcParam::RENS: {
    parameters->setRensMode(mode);
    break;
  }
  case CbcParam::RINS: {
    parameters->setRinsMode(mode);
    break;
  }
  case CbcParam::ROUNDING: {
    parameters->setRoundingMode(mode);
    break;
  }
  case CbcParam::VND: {
    parameters->setVndMode(mode);
    break;
  }
  default:
    break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push a string-valued parameter
*/

int pushCbcSolverStrParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  std::string str = cbcParam.strVal();
  int cbcParamCode = cbcParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (cbcParamCode) {
  case CbcParam::DIRECTORY: {
    char dirSep = CoinFindDirSeparator();
    if (str[str.length() - 1] != dirSep) {
      str += dirSep;
    }
    parameters->setDefaultDirectory(str);
    break;
  }
  default: {
    std::cerr << "pushCbcSolverStr: no equivalent CbcParameters field for "
              << "parameter code `" << cbcParamCode << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  The various parameters to control cut generators can be
  grouped, as they all use the same set of keywords.
*/
int pushCbcSolverCutParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);

  CbcParameters::CGMode mode =
      static_cast<CbcParameters::CGMode>(cbcParam.modeVal());
  int cbcParamCode = cbcParam.paramCode();
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we
     need to do is correct to 0 (no error) if we're successful.
    */
  int retval;

  //if (CoinParamUtils::isInteractive()) {
  //  retval = 1;
  //} else {
    retval = -1;
  //}

  /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
  */
  switch (cbcParamCode) {
  case CbcParam::CLIQUECUTS: {
    parameters->setCliqueMode(mode);
    break;
  }
  case CbcParam::FLOWCUTS: {
    parameters->setFlowMode(mode);
    break;
  }
  case CbcParam::GMICUTS: {
    parameters->setGMIMode(mode);
    break;
  }
  case CbcParam::GOMORYCUTS: {
    parameters->setGomoryMode(mode);
    break;
  }
  case CbcParam::KNAPSACKCUTS: {
    parameters->setKnapsackMode(mode);
    break;
  }
  case CbcParam::LAGOMORYCUTS: {
    parameters->setLaGomoryMode(mode);
    break;
  }
  case CbcParam::LANDPCUTS: {
    parameters->setLandPMode(mode);
    break;
  }
  case CbcParam::LATWOMIRCUTS: {
    parameters->setLaTwomirMode(mode);
    break;
  }
  case CbcParam::MIRCUTS: {
    parameters->setMirMode(mode);
    break;
  }
  case CbcParam::ODDWHEELCUTS: {
    parameters->setOddWheelMode(mode);
    break;
  }
  case CbcParam::PROBINGCUTS: {
    parameters->setProbingMode(mode);
    break;
  }
  case CbcParam::REDSPLITCUTS: {
    parameters->setRedSplitMode(mode);
    break;
  }
  case CbcParam::REDSPLIT2CUTS: {
    parameters->setRedSplit2Mode(mode);
    break;
  }
  case CbcParam::RESIDCAPCUTS: {
    parameters->setResidCapMode(mode);
    break;
  }
  case CbcParam::TWOMIRCUTS: {
    parameters->setTwomirMode(mode);
    break;
  }
  case CbcParam::ZEROHALFCUTS: {
    parameters->setZeroHalfMode(mode);
    break;
  }
  case CbcParam::CUTSTRATEGY: {
    parameters->setCliqueMode(mode);
    parameters->setFlowMode(mode);
    parameters->setGMIMode(mode);
    parameters->setGomoryMode(mode);
    parameters->setKnapsackMode(mode);
    parameters->setLaGomoryMode(mode);
    parameters->setLandPMode(mode);
    parameters->setLaTwomirMode(mode);
    parameters->setMirMode(mode);
    parameters->setOddWheelMode(mode);
    parameters->setProbingMode(mode);
    parameters->setRedSplitMode(mode);
    parameters->setRedSplit2Mode(mode);
    parameters->setRedSplit2Mode(mode);
    parameters->setTwomirMode(mode);
    parameters->setZeroHalfMode(mode);
    break;
  }
  default:
    break;
  }

  return (0);
}

//###########################################################################
//###########################################################################

/*
  This routine imports a new constraint system into the solver.
*/

int doImportParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we
     need to do is correct to 0 (no error) if we're successful.
    */
  int retval;
  //if (CoinParamUtils::isInteractive()) {
  //  retval = 1;
  //} else {
    retval = -1;
  //}
  /*
      Figure out where we're going to acquire this new model. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
  std::string field = cbcParam.strVal();
  std::string fileName;
  if (field == "$") {
     fileName = parameters->getLastMpsIn();
    field = fileName;
  } else if (field == "-") {
    fileName = "stdin";
    field = fileName;
  } else {
    fileName = field;
  }
  /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil
      be the one that actually worked.
    */
  bool canOpen = fileCoinReadable(fileName, parameters->getDefaultDirectory());
  if (canOpen == false) {
    std::cout << "Unable to open file `" << fileName << "', original name '"
              << cbcParam.strVal() << "'." << std::endl;
    return (retval);
  }
  /*
      We can find the file. Record the name. This requires a little finesse:
     what we want is the base file name (and extension(s), if present) but not
     the prefix, unless it's an absolute path.
    */
  if (!fileAbsPath(fileName)) {
    std::string::size_type pos = fileName.rfind(field);
    parameters->setLastMpsIn(fileName.substr(pos));
  } else {
     parameters->setLastMpsIn(fileName);
  }
  /*
      Try to read the file. Standard OSI doesn't support the Clp extensions for
      keepImportNames and allowImportErrors. It should at least support
      keepImportNames. Status will be zero for a successful read.
    */
  OsiSolverInterface *lpSolver = parameters->getModel()->solver();
  int status = lpSolver->readMps(fileName.c_str(), "");
  if (status) {
    std::cout << "There were " << status << " errors on input." << std::endl;
    return (retval);
  }
  /*
      We have a model! Return success.
    */
  parameters->setGoodModel(true);

  return (0);
}

//###########################################################################
//###########################################################################

/*
  This routine imports a debug file into the solver, or arranges for its
  creation. Import works in the standard way, using the file name provided
  with the command.

  As special cases, if the file name is `create' or `createAfterPre', the
  action here sets up to cause a debug file containing the solution to be
  dumped after branch-and-cut is completed.  `createAfterPre' will dump the
  solution before undoing the presolve transforms.  `create' will dump the
  solution after integer presolve is backed out.
*/

int doDebugParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcParameters *parameters = cbcParam.parameters();
  assert(parameters != 0);
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we
     need to do is correct to 0 (no error) if we're successful.
    */
  int retval;
  //if (CoinParamUtils::isInteractive()) {
  //  retval = 1;
  //} else {
    retval = -1;
  //}
  /*
       If the file name is `create' or `createAfterPre', we're just setting up
     to make a debug file the next time we do branch-and-cut.
    */
  std::string field = cbcParam.strVal();
  if (field == "create" || field == "createAfterPre") {
     parameters->setDebugCreate(field);
     return (0);
  }
  /*
      Figure out where we're going to acquire the debug file. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
  std::string fileName;
  if (field == "$") {
     fileName = parameters->getDebugFile();
    field = fileName;
  } else if (field == "-") {
    fileName = "stdin";
    field = fileName;
  } else {
    fileName = field;
  }
  /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil
     be the one that actually worked. No default prefix --- a debug file is
     assumed to always be in the current directory.
    */
  bool canOpen = fileCoinReadable(fileName, parameters->getDefaultDirectory());
  if (canOpen == false) {
    std::cout << "Unable to open file `" << fileName << "', original name '"
              << cbcParam.strVal() << "'." << std::endl;
    return (retval);
  }
  /*
      We can find the file. Record the name. This requires a little finesse:
     what we want is the base file name (and extension(s), if present) but not
     the prefix, unless it's an absolute path.
    */
  if (!fileAbsPath(fileName)) {
    std::string::size_type pos = fileName.rfind(field);
    parameters->setLastMpsIn(fileName.substr(pos));
  } else {
     parameters->setLastMpsIn(fileName);
  }
  /*
      Load the primal variable values into the debug solution vector.
    */
  int intUnused, numCols;
  double dblUnused;
  double *primals;

  bool readOK = readSolution(fileName, intUnused, numCols,
                             dblUnused, 0, 0, &primals, 0);

  if (readOK) {
     parameters->setDebugSol(numCols, primals);
  } else {
    if (primals) {
      delete[] primals;
    }
  }

  return (retval);
}

//###########################################################################
//###########################################################################

void saveSolution(const OsiSolverInterface *osi, std::string fileName)

{
  FILE *fp = fopen(fileName.c_str(), "wb");

  if (fp) {
    int numberRows = osi->getNumRows();
    int numberColumns = osi->getNumCols();
    double objectiveValue = osi->getObjValue();

    fwrite(&numberRows, sizeof(int), 1, fp);
    fwrite(&numberColumns, sizeof(int), 1, fp);
    fwrite(&objectiveValue, sizeof(double), 1, fp);

    const double *primalRowSolution = osi->getRowActivity();
    const double *dualRowSolution = osi->getRowPrice();
    const double *primalColumnSolution = osi->getColSolution();
    const double *dualColumnSolution = osi->getReducedCost();

    fwrite(primalRowSolution, sizeof(double), numberRows, fp);
    fwrite(dualRowSolution, sizeof(double), numberRows, fp);
    fwrite(primalColumnSolution, sizeof(double), numberColumns, fp);
    fwrite(dualColumnSolution, sizeof(double), numberColumns, fp);

    fclose(fp);
  } else {
    std::cout << "saveSolution: Unable to open file `" << fileName << "'."
              << std::endl;
  }

  return;
}

//###########################################################################
//###########################################################################

/*
  Utility routine to read in a solution dump created by saveSolution. Generally
  we don't need all the info in this file, so the routine accepts a bunch of
  reference/pointer paramaters and fills in any that are non-null. It's the
  client's responsibility to dispose of space allocated for solution vectors.
  The parameters fileName, numRows, numCols, and objVal are mandatory. The rest
  can be null.
*/
bool readSolution(std::string fileName, int &numRows, int &numCols,
                  double &objVal, double **rowActivity, double **dualVars,
                  double **primalVars, double **reducedCosts)

{
  FILE *fp = fopen(fileName.c_str(), "rb");
  bool retval = true;

  numRows = -1;
  numCols = -1;
  objVal = 0;
  *rowActivity = 0;
  *dualVars = 0;
  *primalVars = 0;
  *reducedCosts = 0;

  if (fp) {
    fread(&numRows, sizeof(int), 1, fp);
    fread(&numCols, sizeof(int), 1, fp);
    fread(&objVal, sizeof(double), 1, fp);

    if (rowActivity != NULL) {
      *rowActivity = new double[numRows];
      fread(*rowActivity, sizeof(double), numRows, fp);
    } else {
      fseek(fp, numRows * sizeof(double), SEEK_CUR);
    }

    if (dualVars != NULL) {
      *dualVars = new double[numRows];
      fread(*dualVars, sizeof(double), numRows, fp);
    } else {
      fseek(fp, numRows * sizeof(double), SEEK_CUR);
    }

    if (primalVars != NULL) {
      *primalVars = new double[numCols];
      fread(*primalVars, sizeof(double), numCols, fp);
    } else {
      fseek(fp, numCols * sizeof(double), SEEK_CUR);
    }

    if (reducedCosts != NULL) {
      *reducedCosts = new double[numCols];
      fread(*reducedCosts, sizeof(double), numCols, fp);
    } else {
      fseek(fp, numCols * sizeof(double), SEEK_CUR);
    }

    fclose(fp);
  } else {
    std::cout << "readSolution: Unable to open file `" << fileName << "'."
              << std::endl;
    retval = false;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push a double parameter.
*/

int pushCbcModelDblParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcModel *model = cbcParam.model();
  double val = cbcParam.dblVal();
  int cbcParamCode = cbcParam.paramCode();

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcParamCode into the correct key
     for CbcDblParam.
    */
  CbcModel::CbcDblParam key;
  switch (cbcParamCode) {
  case CbcParam::INTEGERTOLERANCE: {
    key = CbcModel::CbcIntegerTolerance;
    break;
  }
  case CbcParam::INFEASIBILITYWEIGHT: {
    key = CbcModel::CbcInfeasibilityWeight;
    break;
  }
  case CbcParam::INCREMENT: {
    key = CbcModel::CbcCutoffIncrement;
    break;
  }
  case CbcParam::ALLOWABLEGAP: {
    key = CbcModel::CbcAllowableGap;
    break;
  }
  case CbcParam::GAPRATIO: {
    key = CbcModel::CbcAllowableFractionGap;
    break;
  }
  case CbcParam::MAXSECONDSNOTIMPROVING: {
    model->setMaxSecondsNotImproving(val);
    break;
  }
  case CbcParam::TIMELIMIT: {
    key = CbcModel::CbcMaximumSeconds;
    break;
  }
  case CbcParam::CUTOFF: {
    key = CbcModel::CbcCurrentCutoff;
    break;
  }
  default: {
    std::cerr << "pushCbcModelDblParam: no equivalent CbcDblParam for "
              << "parameter code `" << cbcParamCode << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  bool result = model->setDblParam(key, val);
  if (result == false) {
    retval = -1;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push an integer parameter.
*/

int pushCbcModelIntParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcModel *model = cbcParam.model();
  int val = cbcParam.intVal();
  int cbcParamCode = cbcParam.paramCode();

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcParamCode into the correct key
     for CbcIntParam, or call the appropriate method directly.
    */
  CbcModel::CbcIntParam key = CbcModel::CbcLastIntParam;
  switch (cbcParamCode) {
  case CbcParam::CUTPASS: {
    model->setMaximumCutPassesAtRoot(val);
    break;
  }
  case CbcParam::LOGLEVEL: {
    CoinMessageHandler *hndl = model->messageHandler();
    assert(hndl != 0);
    hndl->setLogLevel(val);
    break;
  }
  case CbcParam::MAXNODESNOTIMPROVING: {
    model->setMaxNodesNotImproving(val);
    break;
  }
  case CbcParam::MAXSOLS: {
    model->setMaxSolutions(val);
    break;
  }
  case CbcParam::NUMBERBEFORE: {
    model->setNumberBeforeTrust(val);
    break;
  }
  default: {
    std::cerr << "pushCbcModelIntParam: no equivalent CbcIntParam for "
              << "parameter code `" << cbcParamCode << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  if (key != CbcModel::CbcLastIntParam) {
    bool result = model->setIntParam(key, val);
    if (result == false) {
      retval = -1;
    }
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Set CbcModel defaults appropriate for cbc-generic.
*/

void setCbcModelDefaults(CbcModel *model)
{
  model->setIntParam(CbcModel::CbcMaxNumNode, (COIN_INT_MAX / 2));
  model->setIntParam(CbcModel::CbcMaxNumSol, 999999);
  model->setIntParam(CbcModel::CbcFathomDiscipline, 0);

  model->setDblParam(CbcModel::CbcIntegerTolerance, 1.0e-6);
  model->setDblParam(CbcModel::CbcInfeasibilityWeight, 0.0);
  model->setDblParam(CbcModel::CbcCutoffIncrement, 1.0e-5);
  model->setDblParam(CbcModel::CbcAllowableGap, 1.0e-10);
  model->setDblParam(CbcModel::CbcAllowableFractionGap, 0.0);
  // One year is 60x60x24x365 = 31,536,000 seconds.
  model->setDblParam(CbcModel::CbcMaximumSeconds, 3.0e7);
  model->setDblParam(CbcModel::CbcCurrentCutoff, 1.0e100);
  model->setDblParam(CbcModel::CbcOptimizationDirection, 1.0);
  model->setDblParam(CbcModel::CbcCurrentObjectiveValue, 1.0e100);
  model->setDblParam(CbcModel::CbcCurrentMinimizationObjectiveValue, 1.0e100);
  model->setDblParam(CbcModel::CbcStartSeconds, 0.0);

  model->setNumberBeforeTrust(5);
  model->setNumberStrong(5);

  return;
}

} // end namespace CbcParamUtils


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
