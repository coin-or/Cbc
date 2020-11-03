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
  CbcModelParamCode code = cbcParam.paramCode();

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcModelParamCode into the correct key
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
  CbcModelParam &cbcParam = dynamic_cast<CbcModelParam &>(param);
  CbcModel *model = cbcParam.obj();
  int val = cbcParam.intVal();
  CbcModelParamCode code = cbcParam.paramCode();

  assert(model != 0);

  int retval = 0;
  /*
      Translate the parameter code from CbcModelParamCode into the correct key
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

} // end namespace CbcModelParamUtils

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
