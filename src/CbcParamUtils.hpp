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
int doOutputPrecisionParam(CoinParam &param);
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
int pushCbcModelMaxMemoryParam(CoinParam &param);

/** Parse a memory size specification such as "10gb", "500m", "2.5g", a bare
    number of bytes, or one of "unlimited"/"off"/"none" (meaning "no memory
    limit"). Recognized (case-insensitive) unit suffixes are b (bytes,
    default if omitted), k/kb (KiB, 1024 bytes), m/mb (MiB), g/gb (GiB), and
    t/tb (TiB). On success, stores the size in bytes (as a double, since
    byte counts can exceed INT_MAX) in \p bytesOut and returns true;
    returns false (leaving \p bytesOut unchanged) if \p value could not be
    parsed. */
bool parseMemorySize(const std::string &value, double &bytesOut);

void setCbcModelDefaults(CbcModel *model);

} // namespace CbcParamUtils

#endif
