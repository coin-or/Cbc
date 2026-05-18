/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef ClpParamUtils_H
#define ClpParamUtils_H

#include "CoinParam.hpp"
#include "ClpSimplex.hpp"

/*
  Declare the utility functions.
*/

namespace ClpParamUtils {

int doDebugParam(CoinParam &param);
int doExitParam(CoinParam &param);
int doHelpParam(CoinParam &param);
int doImportParam(CoinParam &param);
int doPrintMaskParam(CoinParam &param);
int doNothingParam(CoinParam &param);
int doSolutionParam(CoinParam &param);
int doUnimplementedParam(CoinParam &param);
int doVersionParam(CoinParam &param);

int pushClpDblParam(CoinParam &param);
int pushClpIntParam(CoinParam &param);
int pushClpKwdParam(CoinParam &param);
int pushClpStrParam(CoinParam &param);
int pushClpBoolParam(CoinParam &param);

CLPLIB_EXPORT
void restoreSolution(ClpSimplex *lpSolver, std::string fileName, int mode);
CLPLIB_EXPORT
void saveSolution(const ClpSimplex *lpSolver, std::string fileName);

} // namespace ClpParamUtils

#endif
