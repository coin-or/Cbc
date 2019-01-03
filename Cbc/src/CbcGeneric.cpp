/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include <cassert>
#include <typeinfo>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>

#include "CoinFileIO.hpp"
#include "CoinMpsIO.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

#include "CglCutGenerator.hpp"
#include "CglProbing.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglOddHole.hpp"
#include "CglRedSplit.hpp"
#include "CglTwomir.hpp"

#include "CglPreProcess.hpp"

#include "CbcModel.hpp"
#include "CbcEventHandler.hpp"
#include "CbcTree.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcCompareActual.hpp"

#include "CoinParam.hpp"

#include "CbcGenCtlBlk.hpp"
#include "CbcGenParam.hpp"
#include "CbcGenCbcParam.hpp"
#include "CbcGenOsiParam.hpp"

namespace {

char svnid[] = "$Id$";

}

namespace CbcGenSolvers {
OsiSolverInterface *setupSolvers();
void deleteSolvers();
}

/*
  Unnamed local namespace for cbc-generic support types and functions.
*/

namespace {

/*
  Utility to mark the parameter as having been set by the user. This is a
  bit clumsy --- we need to cast to a derived parameter type to get the
  parameter code. But it'll do 'til I think of a better way.
*/

void markAsSetByUser(CbcGenCtlBlk &ctlBlk, CoinParam *param)

{
  CbcGenParam *genParam = dynamic_cast< CbcGenParam * >(param);
  CbcCbcParam *cbcParam = dynamic_cast< CbcCbcParam * >(param);
  CbcOsiParam *osiParam = dynamic_cast< CbcOsiParam * >(param);
  int code = -1;

  if (genParam != 0) {
    code = genParam->paramCode();
  } else if (cbcParam != 0) {
    code = cbcParam->paramCode();
  } else if (osiParam != 0) {
    code = osiParam->paramCode();
  } else {
    std::cerr
      << "Unrecognised parameter class! Serious internal confusion."
      << std::endl;
  }

  if (code >= 0) {
    ctlBlk.setByUser_[code] = true;
  }

  return;
}

} // end unnamed namespace

int main(int argc, const char *argv[])
{
  /*
      This interior block contains all memory allocation; useful for debugging.
    */
  {
    double time1 = CoinCpuTime();
    double time2;
    /*
          Try and get all the various i/o to come out in order. Synchronise with C
          stdio and make stderr and stdout unbuffered.
        */
    std::ios::sync_with_stdio();
    setbuf(stderr, 0);
    setbuf(stdout, 0);
    /*
          The constructor for ctlBlk establishes the default values for cbc-generic
          parameters. A little more work is required to create the vector of available
          solvers and set the default.
        */
    CbcGenCtlBlk ctlBlk;
    ctlBlk.setMessages();
    ctlBlk.setLogLevel(1);
    OsiSolverInterface *dfltSolver = CbcGenSolvers::setupSolvers();
    ctlBlk.dfltSolver_ = dfltSolver;
    assert(ctlBlk.dfltSolver_);
    /*
          Now we can begin to initialise the parameter vector. Create a vector of the
          proper size, then load up the parameters that are relevant to the main
          program (specifically, values held in ctlBlk and actions evoked from the
          main program).
        */
    int numParams = 0;
    CoinParamVec paramVec;
    paramVec.reserve(CbcOsiParam::CBCOSI_LASTPARAM);
    ctlBlk.paramVec_ = &paramVec;
    ctlBlk.genParams_.first_ = numParams;
    CbcGenParamUtils::addCbcGenParams(numParams, paramVec, &ctlBlk);
    ctlBlk.genParams_.last_ = numParams - 1;
    /*
          Establish a CbcModel object with the default lp solver. Install any defaults
          that are available from ctlBlk.
        */
    CbcModel *model = new CbcModel(*dfltSolver);
    ctlBlk.model_ = model;
    model->messageHandler()->setLogLevel(1);
    model->setNumberStrong(ctlBlk.chooseStrong_.numStrong_);
    model->setNumberBeforeTrust(ctlBlk.chooseStrong_.numBeforeTrust_);
    CbcCbcParamUtils::setCbcModelDefaults(model);
    OsiSolverInterface *osi = model->solver();
    /*
          Set up the remaining classes of parameters, taking defaults from the CbcModel
          and OsiSolverInterface objects we've set up. There are parameters that
          belong to CbcModel (CbcCbcParam) and to the underlying OsiSolverInterface
          (CbcOsiParam).
        */
    ctlBlk.cbcParams_.first_ = numParams;
    CbcCbcParamUtils::addCbcCbcParams(numParams, paramVec, model);
    ctlBlk.cbcParams_.last_ = numParams - 1;
    ctlBlk.osiParams_.first_ = numParams;
    CbcOsiParamUtils::addCbcOsiParams(numParams, paramVec, osi);
    ctlBlk.osiParams_.last_ = numParams - 1;
    /*
          Initialise the vector that tracks parameters that have been changed by user
          command.
        */
    ctlBlk.setByUser_.resize(CbcOsiParam::CBCOSI_LASTPARAM, false);
    /*
          The main command parsing loop. Call getCommand to get the next parameter.
          (The user will be prompted in interactive mode.) If we find something,
          proceed to process it.

          If we don't find anything, behaviour depends on what we've seen so far:

          * An empty command/parameter and no history of previous success gets a
            brief message. If we're in interactive mode, allow the user to try again,
            otherwise quit.

          * An empty command/parameter in interactive mode with some history of
            successful commands is ignored. Iterate and try again.

          * An empty command/parameter when we're not interactive is taken
            as the end of commands. If we have a good model, force branchAndBound.
            (This is one aspect of giving the expected behaviour for
            `cbc-generic [parameters] foo.mps'.)
        */
    bool keepParsing = true;
    bool forceImport = false;
    std::string forceImportFile = "";
    std::string prompt = "cbcGen: ";
    std::string pfx = "";
    while (keepParsing) {
      std::string paramName = CoinParamUtils::getCommand(argc, argv, prompt, &pfx);
      if (paramName.length() == 0) {
        if (ctlBlk.paramsProcessed_ == 0) {
          if (CoinParamUtils::isInteractive()) {
            std::cout
              << "Type `?' or `help' for usage and command keywords."
              << " Type `quit' to quit.";
          } else {
            std::cout
              << "Type `cbc-generic -help' for usage and parameter keywords.";
            keepParsing = false;
          }
          std::cout << std::endl;
        } else if (!CoinParamUtils::isInteractive()) {
          keepParsing = false;
          if (ctlBlk.goodModel_ == true && ctlBlk.bab_.majorStatus_ == CbcGenCtlBlk::BACNotRun) {
            paramName = "branchAndCut";
            pfx = "-";
          }
        }
      }
      if (paramName == "") {
        continue;
      }
      /*
              Do we have a parameter we recognise? In command line mode, if there was no
              prefix (either `-' or `--'), the user didn't intend this as a command
              keyword.
            */
      int matchNdx;
      if (!CoinParamUtils::isCommandLine() || pfx == "-" || pfx == "--") {
        matchNdx = CoinParamUtils::lookupParam(paramName, paramVec);
      } else {
        matchNdx = -3;
      }
      std::cout
        << "Command is `" << paramName
        << "', pfx `" << pfx
        << "', match = " << matchNdx << std::endl;
      /*
              If matchNdx is positive, we have a unique parameter match and we can get on
              with processing. If the return value is negative, and we're not
              interactive, quit. If we're interactive, react as appropriate:
                -1: There was a `?' in the command string. Prompt again.
                -2: No `?', and one or more short matches. Prompt again.
                -3: No `?', but we didn't match anything either. If we're in command line
            	mode, and there was no `-' or `--' prefix, try forcing `import' (but
            	just once, eh). This is the other piece required to get `cbc-generic
            	[parameters] foo.mps' to work as expected.) In interactive mode,
            	we'll require the user to say `import'.  Interactive mode and no
            	history of successful commands gets the help message.
                -4: Configuration error, offer `report to maintainers' message.
            */
      if (matchNdx < 0) {
        if (matchNdx == -3) {
          if (CoinParamUtils::isCommandLine() && pfx == "") {
            if (!forceImport) {
              forceImportFile = paramName;
              paramName = "import";
              matchNdx = CoinParamUtils::lookupParam(paramName, paramVec);
              forceImport = true;
            } else {
              std::cout << "No commands matched `" << paramName << "'."
                        << std::endl;
            }
          } else {
            std::cout << "No commands matched `" << paramName << "'."
                      << std::endl;
            if (ctlBlk.paramsProcessed_ == 0) {
              std::cout
                << "Type `?' or `help' for usage and command keywords."
                << " Type `quit' to quit." << std::endl;
            }
          }
        } else if (matchNdx == -4) {
          std::cout
            << "Please report this error by filing a ticket at "
            << "https://projects.coin-or.org/Cbc/wiki."
            << std::endl;
        }
      }
      if (matchNdx < 0) {
        keepParsing = CoinParamUtils::isInteractive();
        continue;
      }
      CoinParam *param = paramVec[matchNdx];
      ctlBlk.paramsProcessed_++;
      /*
              Depending on the type, we may need a parameter. For keyword parameters, check
              that the keyword is recognised --- setKwdVal will quietly fail on a bad
              keyword.
            */
      CoinParam::CoinParamType type = param->type();
      int valid = 0;
      switch (type) {
      case CoinParam::coinParamAct: {
        break;
      }
      case CoinParam::coinParamInt: {
        int ival = CoinParamUtils::getIntField(argc, argv, &valid);
        if (valid == 0) {
          param->setIntVal(ival);
        }
        break;
      }
      case CoinParam::coinParamDbl: {
        double dval = CoinParamUtils::getDoubleField(argc, argv, &valid);
        if (valid == 0) {
          param->setDblVal(dval);
        }
        break;
      }
      case CoinParam::coinParamStr: {
        if (forceImport) {
          param->setStrVal(forceImportFile);
        } else {
          const std::string tmp = CoinParamUtils::getStringField(argc, argv, &valid);
          if (valid == 0) {
            param->setStrVal(tmp);
          }
        }
        break;
      }
      case CoinParam::coinParamKwd: {
        const std::string tmp = CoinParamUtils::getStringField(argc, argv, &valid);
        if (valid == 0) {
          param->setKwdVal(tmp);
          if (param->kwdVal() != tmp) {
            std::cout
              << "Unrecognised keyword `" << tmp << "' for parameter "
              << param->name() << std::endl;
            param->printKwds();
            std::cout << std::endl;
            valid = 1;
          }
        }
        break;
      }
      default: {
        assert(false);
        break;
      }
      }
      /*
              Deal with missing or incorrect values.

              If valid came back as 2, we're short a parameter. This is interpreted as a
              request to tell the user the current value.  If valid came back as 1, we
              had some sort of parse error. Print an error message.
            */
      if (valid != 0) {
        switch (valid) {
        case 1: {
          std::cout
            << "Could not parse the value given for parameter `"
            << param->name() << "'." << std::endl;
          break;
        }
        case 2: {
          std::cout
            << "Current value of " << param->name() << " parameter is `"
            << *param << "'." << std::endl;
          break;
        }
        default: {
          std::cout
            << "Parse status is " << valid
            << "; this indicates internal confusion." << std::endl
            << "Please report this error by filing a ticket at "
            << "https://projects.coin-or.org/Cbc/wiki."
            << std::endl;
        }
        }
        keepParsing = CoinParamUtils::isInteractive();
        continue;
      }
      /*
              Ok, call the parameter's push function to do the heavy lifting. Push and pull
              functions return 0 for success, 1 for non-fatal error, -1 for fatal error.
            */
      if (param->pushFunc() == 0) {
        std::cout << "Parameter `" << param->name()
                  << "' is not implemented." << std::endl;
      } else {
        int retval = (param->pushFunc())(param);
        markAsSetByUser(ctlBlk, param);
        if (retval < 0) {
          keepParsing = false;
        }
      }
    }
    /*
          End of loop to parse and execute parameter actions. Time to do cleanup.
          The destructor for CbcGenCtlBlk will delete anything with a non-null pointer,
          so we need to be careful that the default solver is deleted only once.
        */
    ctlBlk.dfltSolver_ = 0;
    CbcGenSolvers::deleteSolvers();
    for (int i = 0; i < paramVec.size(); i++) {
      if (paramVec[i] != 0)
        delete paramVec[i];
    }
  }
  /*
      End of memory allocation block. There should be no allocated objects at
      this point.
    */
  return (0);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
