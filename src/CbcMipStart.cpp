#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <OsiSolverInterface.hpp>
#include "CbcMessage.hpp"
#include "CbcHeuristic.hpp"
#include <CbcModel.hpp>
#include "CbcMipStart.hpp"
#include "CbcSOS.hpp"
#include "CbcBoundPropagation.hpp"
#include "CoinTime.hpp"

#define MAX_MIPSTART_WARNINGS_TYPES 10

bool isNumericStr(const char *str)
{
  const size_t l = strlen(str);

  for (size_t i = 0; i < l; ++i)
    if (!(isdigit(str[i]) || (str[i] == '.') || (str[i] == '-') || (str[i] == '+') || (str[i] == 'e')))
      return false;

  return true;
}

std::vector< std::pair< std::string, double > > CbcMipStart::validateMIPStartValues(
  OsiSolverInterface *solver, const std::vector< std::pair< std::string, double > > &colValues,
  CoinMessageHandler *messHandler, CoinMessages *pcoinmsgs,
  bool preProcessedModel)
{
  messHandler->setMessageLimit(CBC_MIPSTART_INVALID_VALUE, MAX_MIPSTART_WARNINGS_TYPES);
  messHandler->setMessageLimit(CBC_MIPSTART_INVALID_COLUMN, MAX_MIPSTART_WARNINGS_TYPES);
  messHandler->setMessageLimit(CBC_MIPSTART_OUT_OF_BOUNDS, MAX_MIPSTART_WARNINGS_TYPES);
  messHandler->setMessageLimit(CBC_MIPSTART_NON_INTEGRAL, MAX_MIPSTART_WARNINGS_TYPES);
  messHandler->setMessageLimit(CBC_MIPSTART_UNKNOWN_COLUMN, MAX_MIPSTART_WARNINGS_TYPES);

  std::vector< std::pair< std::string, double > > validValues;

  CoinMessages &messages = *pcoinmsgs;
  const int numCols = solver->getNumCols();

  std::map< std::string, int > colIndex;
  for (int col = 0; col < numCols; ++col)
    colIndex[solver->getColName(col)] = col;

  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double integerTol = solver->getIntegerTolerance();
  double primalTol = 1e-8;
  solver->getDblParam(OsiPrimalTolerance, primalTol);

  if (!preProcessedModel) {
    // if original model, we would expect all columns to be know
    for (size_t i = 0; i < colValues.size(); ++i) {
      const std::string &colName = colValues[i].first;
      std::map< std::string, int >::const_iterator it = colIndex.find(colName);
      if (it == colIndex.end()) {
        messHandler->message(CBC_MIPSTART_UNKNOWN_COLUMN, messages)
          << colName << colValues[i].second << CoinMessageEol;
      }
    }
  }

  // checking bounds and integrality of each variable and check if value
  // provided in MIP start is Ok
  for (size_t i = 0; i < colValues.size(); ++i) {
    const std::string &colName = colValues[i].first;
    const double value = colValues[i].second;
    std::map< std::string, int >::const_iterator it = colIndex.find(colName);
    if (it == colIndex.end()) {
      // variable does not appears in MIPStart
      continue;
    }

    const int idx = it->second;
    const double lo = lower[idx];
    const double up = upper[idx];
    if (value < lo - primalTol || value > up + primalTol) {
      messHandler->message(CBC_MIPSTART_OUT_OF_BOUNDS, messages)
        << value << colName << lo << up << CoinMessageEol;
      continue;
    }

    bool isInt = solver->isInteger(idx) || solver->isBinary(idx);
    if (isInt) {
      double nearest = floor(value + 0.5);
      double gap = fabs(value - nearest);
      if (gap > integerTol) {
        messHandler->message(CBC_MIPSTART_NON_INTEGRAL, messages)
          << value << colName << gap << integerTol << CoinMessageEol;
        continue;
      }
    }

    validValues.push_back(colValues[i]);
  }

  return validValues;
}

int CbcMipStart::read(OsiSolverInterface *solver, const char *fileName,
  std::vector< std::pair< std::string, double > > &colValues,
  double &solObj, CoinMessageHandler *messHandler, CoinMessages *pcoinmsgs)
{
  CoinMessages &messages = *pcoinmsgs;
#define STR_SIZE 4096
  char printLine[STR_SIZE] = "";
  FILE *f = fopen(fileName, "r");
  if (!f) {
    sprintf(printLine, "Unable to open file %s.", fileName);
    messHandler->message(CBC_MIPSTART_WARNING, messages) << printLine << CoinMessageEol;
    return 1;
  }
  char line[STR_SIZE] = "";

  int nLine = 0;
  // check if psv format!
  int lengthFilename = strlen(fileName);
  // separator
  char separator = ' ';
  if (strstr(fileName, ".psv") == fileName + lengthFilename - 4)
    separator = '|';
  else if (strstr(fileName, ".csv") == fileName + lengthFilename - 4)
    separator = ',';
  if (separator == ' ') {
    /* Two layouts are in circulation and both are accepted, decided per line:

         index name value [reduced-cost]   as written by cbc's own -solution
         name value                        as used by MIPLIB 2017 .sol files

       Deciding per line rather than per file also skips the free-text header
       cbc writes ("Optimal - objective value 3.00000000") without needing to
       recognise it. The one ambiguous shape is a two-field line whose first
       field is numeric: "12 3" could be index 12 with no value, or a column
       literally named "12" holding 3. It is read as the latter, since the
       four-field layout always supplies a value, and a column whose name is a
       number is unusual but legal. */
    while (fgets(line, STR_SIZE, f)) {
      ++nLine;
      char col[4][STR_SIZE] = { "", "", "", "" };
      int nread = sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);
      if (nread <= 0)
        continue;
      const char *name;
      const char *valueStr;
      if (nread >= 3 && strlen(col[0]) && isdigit(col[0][0])) {
        /* index name value ... */
        if (!isNumericStr(col[0])) {
          messHandler->message(CBC_MIPSTART_INVALID_COLUMN, messages) << fileName << nLine << "first" << CoinMessageEol;
          continue;
        }
        name = col[1];
        valueStr = col[2];
      } else if (nread == 2) {
        /* name value */
        name = col[0];
        valueStr = col[1];
      } else {
        /* header or comment - not a variable line */
        continue;
      }
      if (!isNumericStr(valueStr)) {
        messHandler->message(CBC_MIPSTART_INVALID_COLUMN, messages) << fileName << nLine
          << (valueStr == col[2] ? "third" : "second") << CoinMessageEol;
        continue;
      }

      colValues.push_back(std::pair< std::string, double >(std::string(name), atof(valueStr)));
    }
  } else {
    // csv or psv
    int nBad1 = 0;
    int nBad2 = 0;
    while (fgets(line, STR_SIZE, f)) {
      ++nLine;
      // clean line
      // out \n \r and blanks
      int n = strlen(line);
      int nNew = 0;
      for (int i = 0; i < n; i++) {
        char charX = line[i];
        if (charX == ' ') {
          continue;
        } else if (charX == '\n' || charX == '\r') {
          line[nNew] = '\0';
          break;
        } else {
          line[nNew++] = charX;
        }
      }
      char *pipeorcomma = strchr(line, separator);
      if (!pipeorcomma) {
        if (!nBad1) {
          if (nLine > 1) {
            sprintf(printLine, "Reading: %s, line %d (%s) - mipstart file should contain |.", fileName, nLine, line);
            messHandler->message(CBC_MIPSTART_WARNING, messages) << printLine << CoinMessageEol;
          } else {
            // may be OK
            nBad1--;
          }
        }
        nBad1++;
        continue;
      }
      *pipeorcomma = '\0';
      if (!isNumericStr(pipeorcomma + 1)) {
        if (!nBad2) {
          if (nLine > 1) {
            messHandler->message(CBC_MIPSTART_INVALID_COLUMN, messages) << fileName << nLine << "second" << CoinMessageEol;
          } else {
            // may be OK
            nBad2--;
          }
        }
        nBad2++;
        continue;
      }

      double value = atof(pipeorcomma + 1);

      colValues.push_back(std::pair< std::string, double >(std::string(line), value));
    }
    if (nBad1 || nBad2) {
      sprintf(printLine, "Reading: %s, %d errors.", fileName, nBad1 + nBad2);
      messHandler->message(CBC_MIPSTART_WARNING, messages) << printLine << CoinMessageEol;
      return 1;
    }
  }

  const int numCols = solver->getNumCols();
  if (colValues.size()) {
    sprintf(printLine, "MIPStart values read for %d variables.", static_cast< int >(colValues.size()));
    messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
    std::vector< std::pair< std::string, double > > fullValues;
    /* for fast search of column names */
    std::map< std::string, int > colIdx;
    for (int i = 0; i < numCols; i++) {
      fullValues.push_back(std::pair< std::string, double >(solver->getColName(i), 0.0));
      colIdx[solver->getColName(i)] = i;
    }
    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();

    // check bounds and integrality
    std::vector< std::pair< std::string, double > > fixedValues =
      validateMIPStartValues(solver, colValues, messHandler, pcoinmsgs, false);

    int nBadValues = colValues.size() - fixedValues.size();
    for (int i = 0; (i < static_cast< int >(colValues.size())); ++i) {
      std::map< std::string, int >::const_iterator mIt = colIdx.find(colValues[i].first);
      if (mIt != colIdx.end()) {
        const int idx = mIt->second;
        double v = colValues[i].second;
        fullValues[idx].second = v;
      }
    }
    if (nBadValues) {
      sprintf(printLine, "Warning: modifying %d solution values outside bounds",
        nBadValues);
      messHandler->message(CBC_MIPSTART_WARNING, messages) << printLine << CoinMessageEol;
    }
    colValues = fullValues;
  } else {
    sprintf(printLine, "File %s yielded no variable values - expected lines of"
      " \"index name value\" or \"name value\"; MIPStart ignored.", fileName);
    messHandler->message(CBC_MIPSTART_WARNING, messages) << printLine << CoinMessageEol;
    fclose(f);
    return 1;
  }

  fclose(f);
  return 0;
}

int CbcMipStart::computeCompleteSolution(CbcModel *model, OsiSolverInterface *solver,
  const std::vector< std::string > &colNames,
  const std::vector< std::pair< std::string, double > > &colValues,
  double *sol, double &obj, int extraActions, CoinMessageHandler *messHandler, CoinMessages *pmessages,
  int fixMode)
{
  if (!solver->getNumCols())
    return 0;

  const bool fixContinuousToo = (fixMode == FixIntegersAndContinuous);
  const bool assumeUnmentionedZero = (fixMode == FixIntegersAssumeZero);

  CoinMessages &messages = *pmessages;

  int status = 0;
  double compObj = COIN_DBL_MAX;
  bool foundIntegerSol = false;
  OsiSolverInterface *lp = solver->clone();

  std::map< std::string, int > colIdx;
  assert((static_cast< int >(colNames.size())) == lp->getNumCols());
  /* for fast search of column names */
  for (int i = 0; (i < static_cast< int >(colNames.size())); ++i)
    colIdx[colNames[i]] = i;

  char printLine[STR_SIZE];
  int fixed = 0;
  int notFound = 0;
  char colNotFound[256] = "";
  int nContinuousFixed = 0;
  /* The true objective, snapshotted before the recovery LP runs: the cost
     reported for the start has to be its cost in the user's problem, and
     smallBranchAndBound() below is free to reshape the objective of the solver
     it is handed. */
  double *realObj = new double[lp->getNumCols()];
  memcpy(realObj, lp->getObjCoefficients(), sizeof(double) * lp->getNumCols());

  /* There used to be a disabled alternative here that replaced the objective
     with all-ones before re-solving, on the theory that unfixed variables would
     then sink to zero of their own accord and spare us the mini branch and bound
     below. Measured over the 464-instance mip-sanity corpus in all three input
     shapes, it is worse on 36 instances and better on none: it does eliminate
     the only mini-B&B firing (mod011, 16 fractional variables, 0.93 s) but pays
     for it with a solution costing 1.13e+07 instead of -5.46e+07, because the
     LP it solves is no longer the user's problem and lands on whichever feasible
     point is nearest, not a good one. seymour still needed the mini-B&B either
     way. Cheapening the one search we do run, at the price of the answer, is not
     a trade worth keeping, so the alternative is gone rather than dormant. */

  if (assumeUnmentionedZero) {
    // all not mentioned are at zero
    for (int i = 0; (i < lp->getNumCols()); ++i) {
      if (lp->isInteger(i))
        lp->setColBounds(i, 0.0, 0.0);
    }
  }
  if (extraActions) {
    const double *objective = lp->getObjCoefficients();
    const double *lower = lp->getColLower();
    const double *upper = lp->getColUpper();
    for (int i = 0; (i < lp->getNumCols()); ++i) {
      if (lp->isInteger(i)) {
        double objValue = objective[i];
        double lowerValue = lower[i];
        double upperValue = upper[i];
        switch (extraActions) {
        case 1:
          lp->setColBounds(i, lowerValue, lowerValue);
          break;
        case 2:
          lp->setColBounds(i, upperValue, upperValue);
          break;
        case 3:
          lp->setColBounds(i, lowerValue, lowerValue);
          if (objValue < 0.0)
            lp->setColBounds(i, upperValue, upperValue);
          break;
        case 4:
          lp->setColBounds(i, upperValue, upperValue);
          if (objValue > 0.0)
            lp->setColBounds(i, lowerValue, lowerValue);
          break;
        case 5:
          lp->setColBounds(i, lowerValue, lowerValue);
          if (objValue > 0.0)
            lp->setColBounds(i, upperValue, upperValue);
          break;
        case 6:
          lp->setColBounds(i, upperValue, upperValue);
          if (objValue < 0.0)
            lp->setColBounds(i, lowerValue, lowerValue);
          break;
        }
      }
    }
  }
  for (int i = 0; (i < static_cast< int >(colValues.size())); ++i) {
    std::map< std::string, int >::const_iterator mIt = colIdx.find(colValues[i].first);
    if (mIt == colIdx.end()) {
      if (!notFound)
        strcpy(colNotFound, colValues[i].first.c_str());
      notFound++;
    } else {
      const int idx = mIt->second;
      double v = colValues[i].second;
      if (!fixContinuousToo && !lp->isInteger(idx))
        continue;
      if (fabs(v) < 1e-8)
        v = 0.0;
      if (lp->isInteger(idx)) // just to avoid small
        v = floor(v + 0.5); // fractional garbage
      else
        nContinuousFixed++;

      lp->setColBounds(idx, v, v);
      ++fixed;
    }
  }

  if (extraActions)
    fixed = lp->getNumIntegers();
  if (!fixed) {
    // but might be SOS
    if (model) {
      int numberObjects = model->numberObjects();
      for (int i = 0; i < numberObjects; i++) {
        const CbcSOS *object = dynamic_cast< const CbcSOS * >(model->object(i));
        if (object) {
          fixed = 1;
          break; // SOS assume user is expert
        }
      }
    }
    if (!fixed) {
      messHandler->message(CBC_GENERAL, messages)
        << "Warning: MIPstart solution is not valid, column names do not match, ignoring it."
        << CoinMessageEol;
      goto TERMINATE;
    }
  }

  if (notFound >= ((static_cast< double >(colNames.size())) * 0.5)) {
    sprintf(printLine, "Warning: %d column names were not found (e.g. %s) while filling solution.", notFound, colNotFound);
    messHandler->message(CBC_GENERAL, messages)
      << printLine << CoinMessageEol;
  }
  /* Most columns are fixed by now whichever mode is in force, so presolve has a
     lot to remove and the LP left to solve is much smaller. */
  lp->setHintParam(OsiDoPresolveInInitial, true, OsiHintDo);

  // lp->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
  lp->initialSolve();

  if ((lp->isProvenPrimalInfeasible()) || (lp->isProvenDualInfeasible())) {
    /* Say *where* the fixed problem broke, before the recovery below starts
       unfixing things and the evidence is gone. The LP only reports that no
       feasible point exists; bound propagation reasons row by row, so when it
       reaches the same verdict it can name the row (and sometimes the column)
       that closed, which is the one thing a user needs in order to correct a
       start that is genuinely wrong.

       On a clone, and for reporting only. Measured over the mip-sanity corpus
       (464 instances x fix modes "integerZero" and "all"), propagation agreed
       with the LP on all 928 cases: it never contradicted a feasible fixing,
       and in both infeasible ones it proved infeasibility itself and named the
       row. Its *filling* was not worth having on the same evidence -- across
       the 243 cases with columns still free it fixed extra columns in 29, and
       the LP derives those values correctly anyway -- so the tightenings are
       deliberately dropped rather than applied: a start that is merely
       imprecise must keep its chance at the recovery path below.

       Both corpus infeasibilities were exactly that, imprecision: residuals of
       2e-09 (dsbmip) and 1.25e-06 (jssp_la06) against a row, no bound
       violation, and both solve to optimality at a primal tolerance of 1e-5.
       So this branch is reached by rounding noise in a valid solution at least
       as often as by a real error, and the message says "detected at", not
       "caused by". */
    if (messHandler->logLevel() > 0) {
      OsiSolverInterface *bpLp = lp->clone();
      CbcBoundPropagation bp;
      const double bpStart = CoinCpuTime();
      if (!bp.run(bpLp, NULL, 0, CbcBoundPropagation::Fixpoint,
            0, false, 10.0, bpStart)) {
        const int iRow = bp.infeasibleRow();
        const int iCol = bp.infeasibleCol();
        if (iRow >= 0 && iRow < bpLp->getNumRows())
          sprintf(printLine, "Bound propagation confirms the fixed problem is"
            " infeasible, detected at row %d (%s).",
            iRow, bpLp->getRowName(iRow).c_str());
        else if (iCol >= 0 && iCol < bpLp->getNumCols())
          sprintf(printLine, "Bound propagation confirms the fixed problem is"
            " infeasible, detected at column %d (%s).",
            iCol, colNames[iCol].c_str());
        else
          sprintf(printLine, "Bound propagation confirms the fixed problem is"
            " infeasible (no single row or column identified).");
        messHandler->message(CBC_GENERAL, messages)
          << printLine << CoinMessageEol;
      }
      delete bpLp;
    }
    if (nContinuousFixed) {
      messHandler->message(CBC_GENERAL, messages)
        << "Trying just fixing integer variables (and fixingish SOS)." << CoinMessageEol;
      int numberColumns = lp->getNumCols();
      const double *oldLower = solver->getColLower();
      const double *oldUpper = solver->getColUpper();
      double *savedSol = CoinCopyOfArray(lp->getColLower(), numberColumns);
      for (int i = 0; i < numberColumns; ++i) {
        if (!lp->isInteger(i)) {
          lp->setColLower(i, oldLower[i]);
          lp->setColUpper(i, oldUpper[i]);
        }
      }
      // but look at SOS
      if (model) {
        int numberObjects = model->numberObjects();
        for (int i = 0; i < numberObjects; i++) {
          const CbcSOS *object = dynamic_cast< const CbcSOS * >(model->object(i));
          if (object) {
            int n = object->numberMembers();
            const int *members = object->members();
            int sosType = object->sosType();
            if (sosType == 1) {
              // non zero can take any value - others zero
              int iColumn = -1;
              for (int j = 0; j < n; j++) {
                int jColumn = members[j];
                if (savedSol[jColumn])
                  iColumn = jColumn;
              }
              for (int j = 0; j < n; j++) {
                int jColumn = members[j];
                if (jColumn != iColumn) {
                  lp->setColLower(jColumn, 0.0);
                  lp->setColUpper(jColumn, 0.0);
                }
              }
            } else if (sosType == 2) {
              // SOS 2 - make a guess if just one nonzero
              int jA = -1;
              int jB = -1;
              for (int j = 0; j < n; j++) {
                int jColumn = members[j];
                if (savedSol[jColumn]) {
                  if (jA == -1)
                    jA = j;
                  jB = j;
                }
              }
              if (jB > jA + 1) {
                jB = jA + 1;
              } else if (jA == jB) {
                if (jA == n - 1)
                  jA--;
                else
                  jB++;
              }
              for (int j = 0; j < n; j++) {
                if (j != jA && j != jB) {
                  int jColumn = members[j];
                  lp->setColLower(jColumn, 0.0);
                  lp->setColUpper(jColumn, 0.0);
                }
              }
            }
          }
        }
      }

      delete[] savedSol;
      lp->initialSolve();
    } else {
      messHandler->message(CBC_GENERAL, messages)
        << "Fixing only non-zero variables." << CoinMessageEol;
      /* unfix all variables which are zero */
      int notZeroAnymore = 0;
      for (int i = 0; (i < lp->getNumCols()); ++i)
        if (((fabs(lp->getColLower()[i])) <= 1e-8) && (fabs(lp->getColLower()[i] - lp->getColUpper()[i]) <= 1e-8)) {
          const double *oldLower = solver->getColLower();
          const double *oldUpper = solver->getColUpper();
          lp->setColLower(i, oldLower[i]);
          lp->setColUpper(i, oldUpper[i]);
          notZeroAnymore++;
        }
      if (notZeroAnymore)
        lp->initialSolve();
    }
  }

  if (!lp->isProvenOptimal()) {
    messHandler->message(CBC_GENERAL, messages)
      << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
    status = 1;
    goto TERMINATE;
  }
  /* some additional effort is needed to provide an integer solution */
  if (lp->getFractionalIndices().size() > 0) {
    sprintf(printLine, "MIPStart solution provided values for %d of %d integer variables, %d variables are still fractional.", fixed, lp->getNumIntegers(), static_cast< int >(lp->getFractionalIndices().size()));
    messHandler->message(CBC_GENERAL, messages)
      << printLine << CoinMessageEol;
    if (lp->getFractionalIndices().size() < 5) {
      for (int i = 0; i < lp->getFractionalIndices().size(); i++) {
        int iColumn = lp->getFractionalIndices()[i];
        sprintf(printLine, "Variable %d %s has value %g", iColumn,
          colNames[iColumn].c_str(), lp->getColSolution()[iColumn]);
        messHandler->message(CBC_GENERAL, messages)
          << printLine << CoinMessageEol;
      }
    }
    double start = CoinCpuTime();
    CbcSerendipity heuristic(*model);
    /* Anything in (1.0,1000000.0) is read by smallBranchAndBound() as a packed
       "type.fraction" code rather than a plain size ratio: it takes the
       fractional part as the threshold and the tens digit as a node schedule.
       2.0 has fractional part 0.0, so "ratio > fraction" holds for every model
       and the node limit was rewritten to -1 -- no search at all -- discarding
       the 1000 asked for below. On seymour that lost a genuine optimum (423
       reported as unusable, search ended at 430) although a plain
       branch-and-bound on the same LP finds it without branching. 1.0 keeps the
       intended meaning of "any size is small enough" while leaving the node
       limit alone. */
    heuristic.setFractionSmall(1.0);
    heuristic.setFeasibilityPumpOptions(1008013);
    /* 1000 nodes was measured to be too tight: CbcStrategyDefaultSubTree's
       lightweight cut set (essentially just probing, none of the root
       model's Gomory/knapsack/MIR/zero-half cuts) can need well over 1000
       nodes to find *any* feasible completion on an equality-constrained
       fixture no larger than 4 rows x 20 binaries (market split), even
       though the same LP solves in under 100 nodes with the full cut set.
       10000 gives enough headroom to be robust to that without materially
       changing the cost of the common case, where a completion is found in
       a handful of nodes. */
    int returnCode = heuristic.smallBranchAndBound(lp,
      10000, sol,
      compObj,
      model->getCutoff(),
      "ReduceInMIPStart");
    if ((returnCode & 1) != 0) {
      sprintf(printLine, "Mini branch and bound defined values for remaining variables in %.2f seconds.",
        CoinCpuTime() - start);
      messHandler->message(CBC_GENERAL, messages)
        << printLine << CoinMessageEol;
      foundIntegerSol = true;
      lp->getDblParam(OsiObjOffset, obj);
      obj = -obj;
      for (int i = 0; (i < lp->getNumCols()); ++i)
        obj += realObj[i] * sol[i];
    } else {
      messHandler->message(CBC_GENERAL, messages)
        << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
      status = 1;
      goto TERMINATE;
    }
  } else {
    foundIntegerSol = true;

    lp->getDblParam(OsiObjOffset, obj);
    obj = -obj;
    for (int i = 0; (i < lp->getNumCols()); ++i)
      obj += realObj[i] * lp->getColSolution()[i];
    compObj = obj;
    std::copy(lp->getColSolution(), lp->getColSolution() + lp->getNumCols(), sol);
  }

  if (foundIntegerSol) {
    sprintf(printLine, "MIPStart provided solution with cost %g", compObj);
    messHandler->message(CBC_GENERAL, messages)
      << printLine << CoinMessageEol;
    for (int i = 0; (i < lp->getNumCols()); ++i) {
      if (lp->isInteger(i))
        sol[i] = floor(sol[i] + 0.5);
    }
  }

TERMINATE:
  delete[] realObj;
  delete lp;
  return status;
}
#undef STR_SIZE

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
