/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef CbcParamUtils_H
#define CbcParamUtils_H

#include "CoinParam.hpp"

#include "OsiSolverInterface.hpp"

/*
  Declare the utility functions.
*/

namespace CbcParamUtils {

void saveSolution(const OsiSolverInterface *osi, std::string fileName);
bool readSolution(std::string fileName, int &numRows, int &numCols,
                  double &objVal, double **rowActivity, double **dualVars,
                  double **primalVars, double **reducedCosts);

int doBaCParam(CoinParam &param);
int doDebugParam(CoinParam &param);
int doExitParam(CoinParam &param);
int doHelpParam(CoinParam &param);
int doImportParam(CoinParam &param);
int doPrintMaskParam(CoinParam &param);
int doNothingParam(CoinParam &param);
int doSolutionParam(CoinParam &param);
int doUnimplementedParam(CoinParam &param);
int doVersionParam(CoinParam &param);

int pushCbcSolverDblParam(CoinParam &param);
int pushCbcSolverIntParam(CoinParam &param);
int pushCbcSolverKwdParam(CoinParam &param);
int pushCbcSolverStrParam(CoinParam &param);
int pushCbcSolverBoolParam(CoinParam &param);

int pushCbcSolverHeurParam(CoinParam &param);
int pushCbcSolverCutParam(CoinParam &param);

int pushCbcModelDblParam(CoinParam &param);
int pushCbcModelIntParam(CoinParam &param);

void setCbcModelDefaults(CbcModel *model);

} // namespace CbcParamUtils

#endif
