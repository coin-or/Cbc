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
#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcParam.hpp"
#include "CbcParamUtils.hpp"
#include "CbcSettings.hpp"

namespace CbcParamUtils {

//###########################################################################
//###########################################################################

/* Functions to implement  cbc-generic (CbcParam) parameters */

/*
  Maintainer's utility, scan the parameters and report the ones that are
  unimplemented (i.e., have no pushFunc).
*/

int doUnimplementedParam(CoinParam &param)

{
  CbcParam &cbcParam = dynamic_cast<CbcParam &>(param);
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  CoinParamVec &cbcParameters = cbcSettings->getParameters();

  int unimpCnt = 0;
  int maxAcross = 5;
  for (unsigned i = 0; i < cbcParameters.size(); i++) {
    if (cbcParameters[i].pushFunc() == 0) {
      if (unimpCnt % maxAcross == 0) {
        std::cout << std::endl;
      } else {
        std::cout << " ";
      }
      std::cout << cbcParameters[i].name();
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  std::cout << "Cbc version " << cbcSettings->getVersion() << std::endl;

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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int verbose = cbcSettings->getVerbose();
  bool shortHelp = ((verbose & 0x01) ? true : false);
  bool longHelp = ((verbose & 0x02) ? true : false);
  bool hidden = ((verbose & 0x08) ? true : false);

  CoinParamVec &cbcParameters = cbcSettings->getParameters();
  /*
      Tune up the initial cbcSettings. FULLGENERALQUERY will print normally
     hidden params, and a request for long help overrules a request for short
     help.
    */
  if (code == FULLGENERALQUERY) {
    hidden = true;
  }
  if (longHelp) {
    shortHelp = false;
  }

  CoinParamUtils::printGenericHelp();

  std::cout << "\nAvailable commands are:";
  std::string pfx("  ");
  CoinParamUtils::printHelp(cbcParameters, 0, cbcParameters.size() - 1,
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  double val = cbcParam.dblVal();
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
   case ARTIFICIALCOST: {
      cbcSettings->setArtVarMode(CbcSettings::ParamOn, val);
      break;
   }
   case DEXTRA3: {
      cbcSettings->setExtraDbl3(val);
      break;
   }
   case DEXTRA4: {
      cbcSettings->setExtraDbl4(val);
      break;
   }
   case DEXTRA5: {
      cbcSettings->setExtraDbl5(val);
      break;
   }
   case DJFIX: {
      cbcSettings->setDjFixMode(CbcSettings::ParamOn, val);
      break;
   }
   case FAKECUTOFF: {
      cbcSettings->setFeasPumpFakeCutoff(val);
      break;
   }
   case FAKEINCREMENT: {
      cbcSettings->setFeasPumpFakeIncrement(val);
      break;
   }
   case SMALLBAB: {
      cbcSettings->setSmallBab(val);
      break;
   }
   case TIGHTENFACTOR: {
      cbcSettings->setTightenFactor(val);
      break;
   }
     
   default: {
      std::cerr << "pushCbcSolverDbl: no equivalent CbcSettings field for "
                << "parameter code `" << code << "'." << std::endl;
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  int val = cbcParam.intVal();
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
  case BKPIVOTINGSTRATEGY: {
    cbcSettings->setBkPivotStrategy(val);
    break;
  }
  case BKMAXCALLS: {
    cbcSettings->setBkMaxCalls(val);
    break;
  }
  case BKCLQEXTMETHOD: {
    cbcSettings->setBkClqExtMethod(val);
    break;
  }
  case CPP: {
    cbcSettings->setCppMode(val);
    break;
  }
  case CUTDEPTH: {
    cbcSettings->setCutDepth(val);
    break;
  }
  case CUTLENGTH: {
    cbcSettings->setCutLength(val);
    break;
  }
  case CUTPASSINTREE: {
    cbcSettings->setCutPassInTree(val);
    break;
  }
  case DEPTHMINIBAB: {
    cbcSettings->setDepthMiniBaB(val);
    break;
  }
  case DIVEOPT: {
    cbcSettings->setDiveOpt(val);
    break;
  }
  case DIVEOPTSOLVES: {
    cbcSettings->setDiveOptSolves(val);
    break;
  }
  case EXPERIMENT: {
    cbcSettings->setExperimentMode(val);
    break;
  }
  case EXTRA1: {
    cbcSettings->setExtraIntParam1(val);
    break;
  }
  case EXTRA2: {
    cbcSettings->setExtraIntParam2(val);
    break;
  }
  case EXTRA3: {
    cbcSettings->setExtraIntParam3(val);
    break;
  }
  case EXTRA4: {
    cbcSettings->setExtraIntParam4(val);
    break;
  }
  case FPUMPITS: {
    cbcSettings->setFeasPumpIters(val);
    break;
  }
  case FPUMPTUNE: {
    cbcSettings->setFeasPumpTune(val);
    break;
  }
  case FPUMPTUNE2: {
    cbcSettings->setFeasPumpTune2(val);
    break;
  }
  case HEUROPTIONS: {
    cbcSettings->setHeurOptions(val);
    break;
  }
  case LOGLEVEL: {
    cbcSettings->setLogLevel(val);
    break;
  }
  case LPLOGLEVEL: {
    cbcSettings->setLpLogLevel(val);
    break;
  }
  case MAXSAVEDSOLS: {
    cbcSettings->setMaxSavedSols(val);
    break;
  }
  case MAXSLOWCUTS: {
    cbcSettings->setMaxSlowCuts(val);
    break;
  }
  case MOREMOREMIPOPTIONS: {
    cbcSettings->setMoreMoreOptions(val);
    break;
  }
  case MULTIPLEROOTS: {
    cbcSettings->setMultipleRoots(val);
    break;
  }
  case ODDWEXTMETHOD: {
    cbcSettings->setOddWextMethod(val);
    break;
  }
  case OUTPUTFORMAT: {
    cbcSettings->setOutputFormat(val);
    break;
  }
  case PRINTOPTIONS: {
    cbcSettings->setPrintOptions(val);
    break;
  }
  case PROCESSTUNE: {
    cbcSettings->setProcessTune(val);
    break;
  }
  case RANDOMSEED: {
    cbcSettings->setRandomSeed(val);
    break;
  }
  case STRONGSTRATEGY: {
    cbcSettings->setStrongStrategy(val);
    break;
  }
  case TESTOSI: {
    cbcSettings->setTestOsi(val);
    break;
  }
#ifdef CBC_THREADS
  case THREADS: {
    cbcSettings->setThreads(val);
    break;
  }
#endif
  case USERCBC: {
    cbcSettings->setUserCbc(val);
    break;
  }
  case VERBOSE: {
    cbcSettings->setVerbose(val);
    break;
  }
  case VUBTRY: {
    cbcSettings->setVubTry(val);
    break;
  }
  default: {
    std::cerr << "pushCbcSolverInt: no method for storing "
              << "parameter code `" << code << "'." << std::endl;
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  int mode = cbcParam.modeVal();
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
  */

  switch (code) {
  case COMMANDPRINTLEVEL: {
    cbcSettings->setCommandDisplayMode(
                 static_cast<CbcSettings::CommandDisplayMode>(mode));
    break;
  }
  case CLQSTRENGTHENING: {
    cbcSettings->setClqStrMode(static_cast<CbcSettings::ClqStrMode>(mode));
    break;
  }
  case BRANCHPRIORITY: {
    cbcSettings->setBranchingPriorityMode(CbcSettings::BPOff);
    break;
  }
  case CUTOFFCONSTRAINT: {
    cbcSettings->setCutoffMode(static_cast<CbcSettings::CutoffMode>(mode));
    break;
  }
  case INTPRINT: {
    cbcSettings->setIntPrintMode(
        static_cast<CbcSettings::IntPrintMode>(mode));
    break;
  }
  case NODESTRATEGY: {
    cbcSettings->setNodeStrategy(
        static_cast<CbcSettings::NodeStrategy>(mode));
    break;
  }
  case ORBITAL: {
    cbcSettings->setOrbitalStrategy(
        static_cast<CbcSettings::OrbitalStrategy>(mode));
    break;
  }
  case PREPROCESS: {
    cbcSettings->setIPPMode(static_cast<CbcSettings::IPPMode>(mode));
    break;
  }
  case SOSPRIORITIZE: {
    cbcSettings->setSOSStrategy(static_cast<CbcSettings::SOSStrategy>(mode));
    break;
  }
  case STRATEGY: {
    cbcSettings->setStrategyMode(
        static_cast<CbcSettings::StrategyMode>(mode));
    break;
  }
  case TIMEMODE: {
    cbcSettings->setClockType(static_cast<CbcSettings::ClockType>(mode));
    break;
  }
  case USECGRAPH: {
    cbcSettings->setCGraphMode(static_cast<CbcSettings::CGraphMode>(mode));
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  // This is ugly, get keyword and set parameter with it instead.
  CbcSettings::OnOffMode mode =
      static_cast<CbcSettings::OnOffMode>(cbcParam.modeVal());
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int retval = 0;

  switch (code) {
  case CPX: {
    cbcSettings->setCPXMode(mode);
    break;
  }
  case DOHEURISTIC: {
    cbcSettings->setDoHeuristicMode(mode);
    break;
  }
  case ERRORSALLOWED: {
    cbcSettings->setImportErrorsMode(mode);
    break;
  }
  case MESSAGES: {
    cbcSettings->setMessagePrefixMode(mode);
    break;
  }
  case PREPROCNAMES: {
    cbcSettings->setPreProcNamesMode(mode);
    break;
  }
  case SOS: {
    cbcSettings->setSOSMode(mode);
    break;
  }
  case USESOLUTION: {
    cbcSettings->setUseSolutionMode(mode);
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  // This is ugly, get keyword and set parameter with it instead.
  CbcSettings::HeurMode mode =
      static_cast<CbcSettings::HeurMode>(cbcParam.modeVal());
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int retval = 0;

  /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
  */
  switch (code) {
  case HEURISTICSTRATEGY: {
    cbcSettings->setGreedyCoverMode(mode);
    cbcSettings->setGreedyEqualityMode(mode);
    cbcSettings->setCrossoverMode(mode);
    cbcSettings->setDinsMode(mode);
    cbcSettings->setDiveCoefficientMode(mode);
    cbcSettings->setDiveFractionalMode(mode);
    cbcSettings->setDiveGuidedMode(mode);
    cbcSettings->setDiveLineSearchMode(mode);
    cbcSettings->setDivePseudocostMode(mode);
    cbcSettings->setDiveVectorLengthMode(mode);
    cbcSettings->setDWMode(mode);
    cbcSettings->setNaiveHeurMode(mode);
    cbcSettings->setPivotAndFixMode(mode);
#if 0
      cbcSettings->setPivotAndComplementMode(mode);
#endif
    cbcSettings->setProximityMode(mode);
    cbcSettings->setRandRoundMode(mode);
    cbcSettings->setRensMode(mode);
    cbcSettings->setRinsMode(mode);
    cbcSettings->setRoundingMode(mode);
    cbcSettings->setVndMode(mode);
    break;
  }
  case COMBINE: {
    cbcSettings->setCombineMode(mode);
    break;
  }
  case CROSSOVER: {
    cbcSettings->setCrossoverMode(mode);
    break;
  }
  case DINS: {
    cbcSettings->setDinsMode(mode);
    break;
  }
  case DIVINGC: {
    cbcSettings->setDiveCoefficientMode(mode);
    break;
  }
  case DIVINGF: {
    cbcSettings->setDiveFractionalMode(mode);
    break;
  }
  case DIVINGG: {
    cbcSettings->setDiveGuidedMode(mode);
    break;
  }
  case DIVINGL: {
    cbcSettings->setDiveLineSearchMode(mode);
    break;
  }
  case DIVINGP: {
    cbcSettings->setDivePseudocostMode(mode);
    break;
  }
  case DIVINGS: {
    cbcSettings->setDiveRandomMode(mode);
    break;
  }
  case DIVINGV: {
    cbcSettings->setDiveVectorLengthMode(mode);
    break;
  }
  case DW: {
    cbcSettings->setDWMode(mode);
    break;
  }
  case FPUMP: {
    cbcSettings->setFeasPumpMode(mode);
    break;
  }
  case GREEDY: {
    cbcSettings->setGreedyCoverMode(mode);
    cbcSettings->setGreedyEqualityMode(mode);
    break;
  }
  case NAIVE: {
    cbcSettings->setNaiveHeurMode(mode);
    break;
  }
  case PIVOTANDFIX: {
    cbcSettings->setPivotAndFixMode(mode);
    break;
  }
#if 0
   case PIVOTANDCOMPLEMENT: {
      cbcSettings->setPivotAndComplementMode(mode);
      break;
   }
#endif
  case PROXIMITY: {
    cbcSettings->setProximityMode(mode);
    break;
  }
  case RANDROUND: {
    cbcSettings->setRandRoundMode(mode);
    break;
  }
  case RENS: {
    cbcSettings->setRensMode(mode);
    break;
  }
  case RINS: {
    cbcSettings->setRinsMode(mode);
    break;
  }
  case ROUNDING: {
    cbcSettings->setRoundingMode(mode);
    break;
  }
  case VND: {
    cbcSettings->setVndMode(mode);
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  std::string str = cbcParam.strVal();
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
  case DIRECTORY: {
    char dirSep = CoinFindDirSeparator();
    if (str[str.length() - 1] != dirSep) {
      str += dirSep;
    }
    cbcSettings->setDefaultDirectory(str);
    break;
  }
  default: {
    std::cerr << "pushCbcSolverStr: no equivalent CbcSettings field for "
              << "parameter code `" << code << "'." << std::endl;
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);

  CbcSettings::CGMode mode =
      static_cast<CbcSettings::CGMode>(cbcParam.modeVal());
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we
     need to do is correct to 0 (no error) if we're successful.
    */
  int retval;

  if (CoinParamUtils::isInteractive()) {
    retval = 1;
  } else {
    retval = -1;
  }

  /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
  */
  switch (code) {
  case CLIQUECUTS: {
    cbcSettings->setCliqueMode(mode);
    break;
  }
  case FLOWCUTS: {
    cbcSettings->setFlowMode(mode);
    break;
  }
  case GMICUTS: {
    cbcSettings->setGMIMode(mode);
    break;
  }
  case GOMORYCUTS: {
    cbcSettings->setGomoryMode(mode);
    break;
  }
  case KNAPSACKCUTS: {
    cbcSettings->setKnapsackMode(mode);
    break;
  }
  case LAGOMORYCUTS: {
    cbcSettings->setLaGomoryMode(mode);
    break;
  }
  case LANDPCUTS: {
    cbcSettings->setLandPMode(mode);
    break;
  }
  case LATWOMIRCUTS: {
    cbcSettings->setLaTwomirMode(mode);
    break;
  }
  case MIRCUTS: {
    cbcSettings->setMirMode(mode);
    break;
  }
  case ODDWHEELCUTS: {
    cbcSettings->setOddWheelMode(mode);
    break;
  }
  case PROBINGCUTS: {
    cbcSettings->setProbingMode(mode);
    break;
  }
  case REDSPLITCUTS: {
    cbcSettings->setRedSplitMode(mode);
    break;
  }
  case REDSPLIT2CUTS: {
    cbcSettings->setRedSplit2Mode(mode);
    break;
  }
  case RESIDCAPCUTS: {
    cbcSettings->setResidCapMode(mode);
    break;
  }
  case TWOMIRCUTS: {
    cbcSettings->setTwomirMode(mode);
    break;
  }
  case ZEROHALFCUTS: {
    cbcSettings->setZeroHalfMode(mode);
    break;
  }
  case CUTSTRATEGY: {
    cbcSettings->setCliqueMode(mode);
    cbcSettings->setFlowMode(mode);
    cbcSettings->setGMIMode(mode);
    cbcSettings->setGomoryMode(mode);
    cbcSettings->setKnapsackMode(mode);
    cbcSettings->setLaGomoryMode(mode);
    cbcSettings->setLandPMode(mode);
    cbcSettings->setLaTwomirMode(mode);
    cbcSettings->setMirMode(mode);
    cbcSettings->setOddWheelMode(mode);
    cbcSettings->setProbingMode(mode);
    cbcSettings->setRedSplitMode(mode);
    cbcSettings->setRedSplit2Mode(mode);
    cbcSettings->setRedSplit2Mode(mode);
    cbcSettings->setTwomirMode(mode);
    cbcSettings->setZeroHalfMode(mode);
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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we
     need to do is correct to 0 (no error) if we're successful.
    */
  int retval;
  if (CoinParamUtils::isInteractive()) {
    retval = 1;
  } else {
    retval = -1;
  }
  /*
      Figure out where we're going to acquire this new model. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
  std::string field = cbcParam.strVal();
  std::string fileName;
  if (field == "$") {
     fileName = cbcSettings->getLastMpsIn();
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
  bool canOpen = fileCoinReadable(fileName, cbcSettings->getDefaultDirectory());
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
    cbcSettings->setLastMpsIn(fileName.substr(pos));
  } else {
     cbcSettings->setLastMpsIn(fileName);
  }
  /*
      Try to read the file. Standard OSI doesn't support the Clp extensions for
      keepImportNames and allowImportErrors. It should at least support
      keepImportNames. Status will be zero for a successful read.
    */
  OsiSolverInterface *lpSolver = cbcSettings->getModel()->solver();
  int status = lpSolver->readMps(fileName.c_str(), "");
  if (status) {
    std::cout << "There were " << status << " errors on input." << std::endl;
    return (retval);
  }
  /*
      We have a model! Return success.
    */
  cbcSettings->setGoodModel(true);

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
  CbcSettings *cbcSettings = cbcParam.settings();
  assert(cbcSettings != 0);
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we
     need to do is correct to 0 (no error) if we're successful.
    */
  int retval;
  if (CoinParamUtils::isInteractive()) {
    retval = 1;
  } else {
    retval = -1;
  }
  /*
       If the file name is `create' or `createAfterPre', we're just setting up
     to make a debug file the next time we do branch-and-cut.
    */
  std::string field = cbcParam.strVal();
  if (field == "create" || field == "createAfterPre") {
     cbcSettings->setDebugCreate(field);
     return (0);
  }
  /*
      Figure out where we're going to acquire the debug file. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
  std::string fileName;
  if (field == "$") {
     fileName = cbcSettings->getDebugFile();
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
  bool canOpen = fileCoinReadable(fileName, cbcSettings->getDefaultDirectory());
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
    cbcSettings->setLastMpsIn(fileName.substr(pos));
  } else {
     cbcSettings->setLastMpsIn(fileName);
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
     cbcSettings->setDebugSol(numCols, primals);
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
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcParamCode into the correct key
     for CbcDblParam.
    */
  CbcModel::CbcDblParam key;
  switch (code) {
  case INTEGERTOLERANCE: {
    key = CbcModel::CbcIntegerTolerance;
    break;
  }
  case INFEASIBILITYWEIGHT: {
    key = CbcModel::CbcInfeasibilityWeight;
    break;
  }
  case INCREMENT: {
    key = CbcModel::CbcCutoffIncrement;
    break;
  }
  case ALLOWABLEGAP: {
    key = CbcModel::CbcAllowableGap;
    break;
  }
  case GAPRATIO: {
    key = CbcModel::CbcAllowableFractionGap;
    break;
  }
  case MAXSECONDSNOTIMPROVING: {
    model->setMaxSecondsNotImproving(val);
    break;
  }
  case TIMELIMIT: {
    key = CbcModel::CbcMaximumSeconds;
    break;
  }
  case CUTOFF: {
    key = CbcModel::CbcCurrentCutoff;
    break;
  }
  default: {
    std::cerr << "pushCbcModelDblParam: no equivalent CbcDblParam for "
              << "parameter code `" << code << "'." << std::endl;
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
  CbcParamCode code = static_cast<CbcParamCode>(cbcParam.paramCode());

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcParamCode into the correct key
     for CbcIntParam, or call the appropriate method directly.
    */
  CbcModel::CbcIntParam key = CbcModel::CbcLastIntParam;
  switch (code) {
  case CUTPASS: {
    model->setMaximumCutPassesAtRoot(val);
    break;
  }
  case LOGLEVEL: {
    CoinMessageHandler *hndl = model->messageHandler();
    assert(hndl != 0);
    hndl->setLogLevel(val);
    break;
  }
  case MAXNODESNOTIMPROVING: {
    model->setMaxNodesNotImproving(val);
    break;
  }
  case MAXSOLS: {
    model->setMaxSolutions(val);
    break;
  }
  case NUMBERBEFORE: {
    model->setNumberBeforeTrust(val);
    break;
  }
  default: {
    std::cerr << "pushCbcModelIntParam: no equivalent CbcIntParam for "
              << "parameter code `" << code << "'." << std::endl;
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
