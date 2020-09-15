/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
/*
  This file is part of cbc-generic.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <cassert>
#include <string>

#include "CoinFinite.hpp"
#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcModelParam.hpp"
#include "CbcSolverParam.hpp"
#include "CbcSolverSettings.hpp"

namespace {}

/*
  Constructors and destructors

  There's a generic constructor and one for integer, double, keyword, string,
  and action parameters.
*/

/*
  Default constructor.
*/
CbcModelParam::CbcModelParam()
    : CoinParam(), paramCode_(CbcModelParamCode(0)), obj_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for double parameter
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code, std::string name,
                             std::string help, double lower, double upper,
                             double defaultValue, std::string longHelp,
                             CoinDisplayPriority displayPriority)
    : CoinParam(name, help, lower, upper, defaultValue, longHelp,
                displayPriority),
      paramCode_(code), obj_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for integer parameter
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code, std::string name,
                             std::string help, int lower, int upper,
                             int defaultValue, std::string longHelp,
                             CoinDisplayPriority displayPriority)
    : CoinParam(name, help, lower, upper, defaultValue, longHelp,
                displayPriority),
      paramCode_(code), obj_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for keyword parameter.
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code, std::string name,
                             std::string help, std::string defaultKwd,
                             int defaultMode, std::string longHelp,
                             CoinDisplayPriority displayPriority)
    : CoinParam(name, help, defaultKwd, defaultMode, longHelp, displayPriority),
      paramCode_(code), obj_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for string parameter.
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code, std::string name,
                             std::string help, std::string defaultValue,
                             std::string longHelp,
                             CoinDisplayPriority displayPriority)
    : CoinParam(name, help, defaultValue, longHelp, displayPriority),
      paramCode_(code), obj_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for action parameter.
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code, std::string name,
                             std::string help, std::string longHelp,
                             CoinDisplayPriority displayPriority)
    : CoinParam(name, help, longHelp, displayPriority), paramCode_(code),
      obj_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Copy constructor.
*/
CbcModelParam::CbcModelParam(const CbcModelParam &orig)
    : CoinParam(orig), paramCode_(orig.paramCode_), obj_(orig.obj_) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Clone
*/

CbcModelParam *CbcModelParam::clone() { return (new CbcModelParam(*this)); }

CbcModelParam &CbcModelParam::operator=(const CbcModelParam &rhs) {
  if (this != &rhs) {
    CoinParam::operator=(rhs);

    paramCode_ = rhs.paramCode_;
    obj_ = rhs.obj_;
  }

  return *this;
}

//###########################################################################
//###########################################################################

/*
  Destructor
*/
CbcModelParam::~CbcModelParam() { /* Nothing more to do */
}

//###########################################################################
//###########################################################################

namespace CbcModelParamUtils {

/* Function to set up cbc (CbcModel) parameters.  */

void addCbcModelParams(CoinParamVec &parameters, CbcModel *model)

{
  for (int code = CbcModelParam::CBCMODEL_FIRSTPARAM + 1;
       code < CbcModelParam::CBCMODEL_LASTPARAM; code++) {
    CoinParam &param = parameters[code];
    static_cast<CbcModelParam &>(param).setObj(model);
  }

  parameters[CbcModelParam::ALLOWABLEGAP].setup(
      "allow!ableGap",
      "Stop when gap between best possible and incumbent is less than this",
      0.0, 1.0e20, 0.0,
      "If the gap between best solution and best possible solution is less "
      "than this then the search will be terminated. Also see ratioGap.");
  parameters[CbcModelParam::ALLOWABLEGAP].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::CUTOFF].setup(
      "cuto!ff", "All solutions must be better than this", -1.0e60, 1.0e60,
      1.0e50,
      "All solutions must be better than this value (in a minimization sense). "
      " This is also set by cbc whenever it obtains a solution and is set to "
      "the value of the objective for the solution minus the cutoff "
      "increment.");
  parameters[CbcModelParam::CUTOFF].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::DIRECTION].setup(
      "direction", "Minimize or maximize", "min!imize",
      CbcModelParam::OptDirMinimize,
      "The default is minimize - use 'direction maximize' for "
      "maximization.\nYou can also use the parameters 'maximize' or "
      "'minimize'.");
  parameters[CbcModelParam::DIRECTION].appendKwd("max!imize",
                                                 CbcModelParam::OptDirMaximize);
  parameters[CbcModelParam::DIRECTION].appendKwd("zero",
                                                 CbcModelParam::OptDirZero);

  parameters[CbcModelParam::INCREMENT].setup(
      "inc!rement",
      "A new solution must be at least this much better than the incumbent",
      -1.0e20, 1.0e20, model->getDblParam(CbcModel::CbcCutoffIncrement),
      "Whenever a solution is found the bound on future solutions is set to "
      "the objective of the solution (in a minimization sense) plus the "
      "specified increment.  If this option is not specified, the code will "
      "try and work out an increment.  E.g., if all objective coefficients are "
      "multiples of 0.01 and only integer variables have entries in objective "
      "then the increment can be set to 0.01.  Be careful if you set this "
      "negative!");
  parameters[CbcModelParam::INCREMENT].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::INFEASIBILITYWEIGHT].setup(
      "inf!easibilityWeight",
      "Each integer infeasibility is expected to cost this much", 0.0, 1.0e20,
      model->getDblParam(CbcModel::CbcInfeasibilityWeight),
      "A primitive way of deciding which node to explore next.  Satisfying "
      "each integer infeasibility is expected to cost this much.");
  parameters[CbcModelParam::INFEASIBILITYWEIGHT].setPushFunc(
      pushCbcModelDblParam);

  parameters[CbcModelParam::INTEGERTOLERANCE].setup(
      "integerT!olerance",
      "For an optimal solution, no integer variable may be farther than this "
      "from an integer value",
      1.0e-20, 0.5, model->getDblParam(CbcModel::CbcIntegerTolerance),
      "When checking a solution for feasibility, if the difference between the "
      "value of a variable and the nearest integer is less than the integer "
      "tolerance, the value is considered to be integral. Beware of setting "
      "this smaller than the primal tolerance.");
  parameters[CbcModelParam::INTEGERTOLERANCE].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::LOGLEVEL].setup(
      "bclog!Level", "Level of detail in Coin branch and Cut output", -1, 63,
      model->messageHandler()->logLevel(),
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");
  parameters[CbcModelParam::LOGLEVEL].setPushFunc(pushCbcModelIntParam);

  parameters[CbcModelParam::MAXIMIZE].setup(
      "max!imize", "Set optimization direction to maximize",
      "The default is minimize - use 'maximize' for maximization.\n A synonym "
      "for 'direction maximize'.",
      CoinParam::displayPriorityHigh);

  parameters[CbcModelParam::MAXNODES].setup(
      "maxN!odes", "Maximum number of nodes to evaluate", 1, 2147483647, 1,
      "This is a repeatable way to limit search.  Normally using time is "
      "easier but then the results may not be repeatable.");
  parameters[CbcModelParam::MAXNODES].setPushFunc(pushCbcModelIntParam);

  parameters[CbcModelParam::MAXNODESNOTIMPROVING].setup(
      "maxNNI!FS",
      "Maximum number of nodes to be processed without improving the incumbent "
      "solution.",
      -1, COIN_INT_MAX, -1,
      "This criterion specifies that when a feasible solution is available, "
      "the search should continue only if better feasible solutions were "
      "produced in the last nodes.");
  parameters[CbcModelParam::MAXNODESNOTIMPROVING].setPushFunc(
      pushCbcModelIntParam);

  parameters[CbcModelParam::MAXSECONDSNOTIMPROVING].setup(
      "secni!fs", "maximum seconds without improving the incumbent solution",
      -1.0, COIN_DBL_MAX, -1.0,
      "With this stopping criterion, after a feasible solution is found, the "
      "search should continue only if the incumbent solution was updated "
      "recently, the tolerance is specified here. A discussion on why this "
      "criterion can be useful is included here: "
      "https://yetanothermathprogrammingconsultant.blogspot.com/2019/11/"
      "mip-solver-stopping-criteria.html .");
  parameters[CbcModelParam::MAXSECONDSNOTIMPROVING].setPushFunc(
      pushCbcModelDblParam);

  parameters[CbcModelParam::MAXSOLS].setup(
      "maxSo!lutions", "Maximum number of feasible solutions to get", 1,
      COIN_INT_MAX, COIN_INT_MAX,
      "You may want to stop after (say) two solutions or an hour. This is "
      "checked every node in tree, so it is possible to get more solutions "
      "from heuristics.");
  parameters[CbcModelParam::MAXSOLS].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::MINIMIZE].setup(
      "min!imize", "Set optimization direction to minimize",
      "The default is minimize - use 'maximize' for maximization.\nThis should "
      "only be necessary if you have previously set maximization. A synonym "
      "for 'direction minimize'.");

  parameters[CbcModelParam::MIPOPTIONS].setup(
      "mipO!ptions", "Dubious options for mip", 0, COIN_INT_MAX, 0, "",
      CoinParam::displayPriorityNone);

  parameters[CbcModelParam::MOREMIPOPTIONS].setup(
      "more!MipOptions", "More dubious options for mip", -1, COIN_INT_MAX, 0,
      "", CoinParam::displayPriorityNone);

#if 0
  parameters[CbcModelParam::NUMBERMINI].setup("miniT!ree",
                                               "Size of fast mini tree", 0, COIN_INT_MAX, 0,
                                               "The idea is that I can do a small tree fast. This is a first try and will hopefully become more sophisticated.",
                                               CoinParam::displayPriorityNone);
#endif

  parameters[CbcModelParam::NUMBERANALYZE].setup(
      "numberA!nalyze", "Number of analysis iterations", -COIN_INT_MAX,
      COIN_INT_MAX, 0,
      "This says how many iterations to spend at the root node analyzing the "
      "problem.  This is a first try and will hopefully become more "
      "sophisticated.",
      CoinParam::displayPriorityNone);

  parameters[CbcModelParam::CUTPASS].setup(
      "passC!uts", "Number of cut passes at root node", -999999, 999999,
      model->getMaximumCutPassesAtRoot(),
      "The default is 100 passes if less than 500 columns, 100 passes (but "
      "stop if the drop is small) if less than 5000 columns, 20 otherwise.");

  parameters[CbcModelParam::GAPRATIO].setup(
      "ratio!Gap",
      "Stop when the gap between the best possible solution and the incumbent "
      "is less than this fraction of the larger of the two",
      0.0, 1.0e20, model->getDblParam(CbcModel::CbcAllowableFractionGap),
      "If the gap between the best solution and the best possible solution is "
      "less than this fraction of the objective value at the root node then "
      "the search will terminate.  See 'allowableGap' for a way of using "
      "absolute value rather than fraction.");
  parameters[CbcModelParam::GAPRATIO].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::TIMELIMIT].setup(
      "sec!onds", "Maximum seconds for branch and cut", -1.0, 1.0e12, -1.0,
      "After this many seconds the program will act as if maximum nodes had "
      "been reached.");
  parameters[CbcModelParam::TIMELIMIT].setPushFunc(pushCbcModelDblParam);

  parameters[CbcModelParam::STRONGBRANCHING].setup(
      "strong!Branching", "Number of variables to look at in strong branching",
      0, 999999, model->numberStrong(),
      "In order to decide which variable to branch on, the code will choose up "
      "to this number of unsatisfied variables and try mini up and down "
      "branches.  The most effective one is chosen. If a variable is branched "
      "on many times then the previous average up and down costs may be used - "
      "see number before trust.");

  parameters[CbcModelParam::NUMBERBEFORE].setup(
      "trust!PseudoCosts", "Number of branches before we trust pseudocosts", -1,
      2000000, model->numberBeforeTrust(),
      "Using strong branching computes pseudo-costs.  After this many times "
      "for a variable we just trust the pseudo costs and do not do any more "
      "strong branching.");
  parameters[CbcModelParam::NUMBERBEFORE].setPushFunc(pushCbcModelIntParam);
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

//###########################################################################
//###########################################################################

/*
  Function to push a double parameter.
*/

int pushCbcModelDblParam(CoinParam &param)

{
  CbcModelParam &cbcParam = dynamic_cast<CbcModelParam &>(param);
  CbcModel *model = cbcParam.obj();
  double val = cbcParam.dblVal();
  CbcModelParam::CbcModelParamCode code = cbcParam.paramCode();

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcModelParamCode into the correct key
     for CbcDblParam.
    */
  CbcModel::CbcDblParam key;
  switch (code) {
  case CbcModelParam::INTEGERTOLERANCE: {
    key = CbcModel::CbcIntegerTolerance;
    break;
  }
  case CbcModelParam::INFEASIBILITYWEIGHT: {
    key = CbcModel::CbcInfeasibilityWeight;
    break;
  }
  case CbcModelParam::INCREMENT: {
    key = CbcModel::CbcCutoffIncrement;
    break;
  }
  case CbcModelParam::ALLOWABLEGAP: {
    key = CbcModel::CbcAllowableGap;
    break;
  }
  case CbcModelParam::GAPRATIO: {
    key = CbcModel::CbcAllowableFractionGap;
    break;
  }
  case CbcModelParam::MAXSECONDSNOTIMPROVING: {
    model->setMaxSecondsNotImproving(val);
    break;
  }
  case CbcModelParam::TIMELIMIT: {
    key = CbcModel::CbcMaximumSeconds;
    break;
  }
  case CbcModelParam::CUTOFF: {
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
  CbcModelParam &cbcParam = dynamic_cast<CbcModelParam &>(param);
  CbcModel *model = cbcParam.obj();
  int val = cbcParam.intVal();
  CbcModelParam::CbcModelParamCode code = cbcParam.paramCode();

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcModelParamCode into the correct key
     for CbcIntParam, or call the appropriate method directly.
    */
  CbcModel::CbcIntParam key = CbcModel::CbcLastIntParam;
  switch (code) {
  case CbcModelParam::CUTPASS: {
    model->setMaximumCutPassesAtRoot(val);
    break;
  }
  case CbcModelParam::LOGLEVEL: {
    CoinMessageHandler *hndl = model->messageHandler();
    assert(hndl != 0);
    hndl->setLogLevel(val);
    break;
  }
  case CbcModelParam::MAXNODESNOTIMPROVING: {
    model->setMaxNodesNotImproving(val);
    break;
  }
  case CbcModelParam::MAXSOLS: {
    model->setMaxSolutions(val);
    break;
  }
  case CbcModelParam::NUMBERBEFORE: {
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

} // end namespace CbcModelParamUtils

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
