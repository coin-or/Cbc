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

#include "ClpParameters.hpp"
#include "ClpParamUtils.hpp"
#include "ClpSimplex.hpp"

namespace ClpParamUtils {

//###########################################################################
//###########################################################################

/* Functions to perform actions related to setting parameters */

/*
  Maintainer's utility, scan the parameters and report the ones that are
  unimplemented (i.e., have no pushFunc).
*/

int doUnimplementedParam(CoinParam &param)

{
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters parameters = *clpParam.parameters();

  int unimpCnt = 0;
  int maxAcross = 5;
  for (int code = ClpParam::FIRSTPARAM + 1;
       code < ClpParam::LASTPARAM; code++) {
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
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
  assert(parameters != 0);

  std::cout << "Clp version " << parameters->getVersion() << std::endl;

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
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
  assert(parameters != 0);

  int clpParamCode = clpParam.paramCode();

  int verbose;
  parameters->getParamVal(ClpParam::VERBOSE, verbose);
  bool shortHelp = ((verbose & 0x01) ? true : false);
  bool longHelp = ((verbose & 0x02) ? true : false);
  bool hidden = ((verbose & 0x08) ? true : false);

  CoinParamVec &paramVec = parameters->paramVec();
  /*
     Tune up the initial clp parameters. FULLGENERALQUERY will print normally
     hidden params, and a request for long help overrules a request for short
     help.
    */
  if (clpParamCode == ClpParam::FULLGENERALQUERY) {
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
  This routine imports a new constraint system into the solver.
*/

int doImportParam(CoinParam &param)

{
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
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

  std::string fileName = clpParam.strVal();

  /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil
      be the one that actually worked.
    */
  bool canOpen = fileCoinReadable(fileName, parameters->getDefaultDirectory());
  if (canOpen == false) {
    std::cout << "Unable to open file `" << fileName << "', original name '"
              << clpParam.strVal() << "'." << std::endl;
    return (retval);
  }

  /*
      Try to read the file. Standard OSI doesn't support the Clp extensions for
      keepImportNames and allowImportErrors. It should at least support
      keepImportNames. Status will be zero for a successful read.
    */
  int status = parameters->getModel()->readMps(fileName.c_str(), "");
  if (status) {
    std::cout << "There were " << status << " errors on input." << std::endl;
    return (retval);
  }
  /*
      We have a model! Return success.
    */
  // TODO What should we do here?
  // parameters->setGoodModel(true);

  return (0);
}

//###########################################################################
//###########################################################################

/*
  Function to push a keyword-valued parameter. This is the catch-all function
  for keyword parameters that don't belong to any other useful grouping.
*/

int pushClpKwdParam(CoinParam &param) {
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
  assert(parameters != 0);

  int mode = clpParam.modeVal();
  int clpParamCode = clpParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
  */

  switch (clpParamCode) {
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

int pushClpBoolParam(CoinParam &param) {
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
  assert(parameters != 0);

  // This is ugly, get keyword and set parameter with it instead.
  ClpParameters::OnOffMode mode =
     static_cast<ClpParameters::OnOffMode>(clpParam.modeVal());
  int clpParamCode = clpParam.paramCode();

  int retval = 0;

  switch (clpParamCode) {
  default:
    break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/*
  Function to push a double-valued parameter.
*/

int pushClpDblParam(CoinParam &param)

{
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
  assert(parameters != 0);

  double val = clpParam.dblVal();
  int clpParamCode = clpParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (clpParamCode) {
   default:
     break;
  }

  return (retval);

}

//###########################################################################
//###########################################################################

/*
  Function to push an integer-valued parameter.
*/

int pushClpIntParam(CoinParam &param)
{
  ClpParam &clpParam = dynamic_cast<ClpParam &>(param);
  ClpParameters *parameters = clpParam.parameters();
  assert(parameters != 0);

  int val = clpParam.intVal();
  int clpParamCode = clpParam.paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (clpParamCode) {
  default:
    break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################

/* Restore a solution from file.
   mode 0 normal, 1 swap rows and columns and primal and dual
   if 2 set then also change signs
*/
void restoreSolution(ClpSimplex *lpSolver, std::string fileName, int mode)
{
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (fp) {
    int numberRows = lpSolver->numberRows();
    int numberColumns = lpSolver->numberColumns();
    int numberRowsFile;
    int numberColumnsFile;
    double objectiveValue;
    size_t nRead;
    nRead = fread(&numberRowsFile, sizeof(int), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    nRead = fread(&numberColumnsFile, sizeof(int), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    nRead = fread(&objectiveValue, sizeof(double), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    double *dualRowSolution = lpSolver->dualRowSolution();
    double *primalRowSolution = lpSolver->primalRowSolution();
    double *dualColumnSolution = lpSolver->dualColumnSolution();
    double *primalColumnSolution = lpSolver->primalColumnSolution();
    if (mode) {
      // swap
      int k = numberRows;
      numberRows = numberColumns;
      numberColumns = k;
      double *temp;
      temp = dualRowSolution;
      dualRowSolution = primalColumnSolution;
      primalColumnSolution = temp;
      temp = dualColumnSolution;
      dualColumnSolution = primalRowSolution;
      primalRowSolution = temp;
    }
    if (numberRows > numberRowsFile || numberColumns > numberColumnsFile) {
      std::cout << "Mismatch on rows and/or columns - giving up" << std::endl;
    } else {
      lpSolver->setObjectiveValue(objectiveValue);
      if (numberRows == numberRowsFile && numberColumns == numberColumnsFile) {
        nRead = fread(primalRowSolution, sizeof(double), numberRows, fp);
        if (nRead != static_cast< size_t >(numberRows))
          throw("Error in fread");
        nRead = fread(dualRowSolution, sizeof(double), numberRows, fp);
        if (nRead != static_cast< size_t >(numberRows))
          throw("Error in fread");
        nRead = fread(primalColumnSolution, sizeof(double), numberColumns, fp);
        if (nRead != static_cast< size_t >(numberColumns))
          throw("Error in fread");
        nRead = fread(dualColumnSolution, sizeof(double), numberColumns, fp);
        if (nRead != static_cast< size_t >(numberColumns))
          throw("Error in fread");
      } else {
        std::cout << "Mismatch on rows and/or columns - truncating" << std::endl;
        double *temp = new double[std::max(numberRowsFile, numberColumnsFile)];
        nRead = fread(temp, sizeof(double), numberRowsFile, fp);
        if (nRead != static_cast< size_t >(numberRowsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberRows, primalRowSolution);
        nRead = fread(temp, sizeof(double), numberRowsFile, fp);
        if (nRead != static_cast< size_t >(numberRowsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberRows, dualRowSolution);
        nRead = fread(temp, sizeof(double), numberColumnsFile, fp);
        if (nRead != static_cast< size_t >(numberColumnsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberColumns, primalColumnSolution);
        nRead = fread(temp, sizeof(double), numberColumnsFile, fp);
        if (nRead != static_cast< size_t >(numberColumnsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberColumns, dualColumnSolution);
        delete[] temp;
      }
      if (mode == 3) {
        int i;
        for (i = 0; i < numberRows; i++) {
          primalRowSolution[i] = -primalRowSolution[i];
          dualRowSolution[i] = -dualRowSolution[i];
        }
        for (i = 0; i < numberColumns; i++) {
          primalColumnSolution[i] = -primalColumnSolution[i];
          dualColumnSolution[i] = -dualColumnSolution[i];
        }
      }
    }
    fclose(fp);
  } else {
    std::cout << "Unable to open file " << fileName << std::endl;
  }
}

// Dump a solution to file
void saveSolution(const ClpSimplex *lpSolver, std::string fileName)
{
  if (strstr(fileName.c_str(), "_fix_read_")) {
    FILE *fp = fopen(fileName.c_str(), "rb");
    if (fp) {
      ClpSimplex *solver = const_cast< ClpSimplex * >(lpSolver);
      restoreSolution(solver, fileName, 0);
      // fix all
      int logLevel = solver->logLevel();
      int iColumn;
      int numberColumns = solver->numberColumns();
      double *primalColumnSolution = solver->primalColumnSolution();
      double *columnLower = solver->columnLower();
      double *columnUpper = solver->columnUpper();
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = primalColumnSolution[iColumn];
        if (value > columnUpper[iColumn]) {
          if (value > columnUpper[iColumn] + 1.0e-6 && logLevel > 1)
            printf("%d value of %g - bounds %g %g\n",
              iColumn, value, columnLower[iColumn], columnUpper[iColumn]);
          value = columnUpper[iColumn];
        } else if (value < columnLower[iColumn]) {
          if (value < columnLower[iColumn] - 1.0e-6 && logLevel > 1)
            printf("%d value of %g - bounds %g %g\n",
              iColumn, value, columnLower[iColumn], columnUpper[iColumn]);
          value = columnLower[iColumn];
        }
        columnLower[iColumn] = value;
        columnUpper[iColumn] = value;
      }
      return;
    }
  }
  FILE *fp = fopen(fileName.c_str(), "wb");
  if (fp) {
    int numberRows = lpSolver->numberRows();
    int numberColumns = lpSolver->numberColumns();
    double objectiveValue = lpSolver->objectiveValue();
    size_t nWrite;
    nWrite = fwrite(&numberRows, sizeof(int), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    nWrite = fwrite(&numberColumns, sizeof(int), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    nWrite = fwrite(&objectiveValue, sizeof(double), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    double *dualRowSolution = lpSolver->dualRowSolution();
    double *primalRowSolution = lpSolver->primalRowSolution();
    nWrite = fwrite(primalRowSolution, sizeof(double), numberRows, fp);
    if (nWrite != static_cast< size_t >(numberRows))
      throw("Error in fwrite");
    nWrite = fwrite(dualRowSolution, sizeof(double), numberRows, fp);
    if (nWrite != static_cast< size_t >(numberRows))
      throw("Error in fwrite");
    double *dualColumnSolution = lpSolver->dualColumnSolution();
    double *primalColumnSolution = lpSolver->primalColumnSolution();
    nWrite = fwrite(primalColumnSolution, sizeof(double), numberColumns, fp);
    if (nWrite != static_cast< size_t >(numberColumns))
      throw("Error in fwrite");
    nWrite = fwrite(dualColumnSolution, sizeof(double), numberColumns, fp);
    if (nWrite != static_cast< size_t >(numberColumns))
      throw("Error in fwrite");
    fclose(fp);
  } else {
    std::cout << "Unable to open file " << fileName << std::endl;
  }
}

}
