// Copyright (C) 2026 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcCoefficientStrengthening.hpp"

#include "CoinMessageHandler.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"

#include <cmath>
#include <vector>

// Reference solution loaded by -debugCuts; set before applyLpMethod() so that
// pre-root-LP strengthening can check itself before the OsiRowCutDebugger is
// active. Same globals CbcBoundPropagation.cpp falls back on -- see the
// comment there for why the debugger itself is not yet attached this early
// when driven from the "-debugCuts" CLI switch (it is already active by this
// point when driven through Cbc_activateRowCutDebugger()/mip-debug-cuts).
extern double *debugSolution;
extern int debugNumberColumns;

namespace {

/// Anything at or beyond this magnitude is treated as infinite, matching the
/// convention used throughout Osi/Cgl (see e.g. CglPreProcess).
const double COEFSTR_INFINITY = 1.0e20;

/// Feasibility tolerance. Same role as papilo's Num::feastol: a slack below
/// this is not considered a real slack, and a value this close to its ceiling
/// is snapped up to it.
const double COEFSTR_FEASTOL = 1.0e-6;

/// A column whose bounds are this close together is fixed; tightening its
/// coefficient is exact but pointless, so it is skipped (as papilo does).
const double COEFSTR_EPSILON = 1.0e-9;

/// One pending right-hand-side change, applied only once the matrix
/// replacement has been verified to have taken effect.
struct RowBoundChange {
  int row;
  double lower;
  double upper;
};

} // end anonymous namespace

CbcCoefficientStrengthening::CbcCoefficientStrengthening()
  : nCoefficients_(0)
  , nRows_(0)
  , timeUsed_(0.0)
{
}

bool CbcCoefficientStrengthening::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel)
{
  nCoefficients_ = 0;
  nRows_ = 0;
  timeUsed_ = 0.0;

  if (!solver || solver->getNumRows() == 0 || solver->getNumCols() == 0)
    return false;

  const double t0 = CoinWallclockTime();

  const int numberRows = solver->getNumRows();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();

  // Scan the solver's own row copy read-only, and only take a mutable copy of
  // the matrix once a change is actually known to be needed. On this corpus
  // the pass finds nothing on 387 of 471 instances, and copying the matrix
  // eagerly made those instances 57% of its total cost -- dt_optimization
  // (1.0M nonzeros, no change) spent 5.7 ms of its 5.9 ms purely on a copy
  // that was then discarded. The scan itself is cheap: it reads bounds and
  // coefficients and writes nothing.
  const CoinPackedMatrix *rowCopy = solver->getMatrixByRow();
  const double *element = rowCopy->getElements();
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();

  // Deferred mutable copy plus the writes destined for it. A change is
  // recorded here during the scan and replayed into the copy afterwards,
  // which keeps the copy off the no-change path entirely.
  //
  // Note this must be a copy of the *pre-change* matrix, so it has to be
  // taken before any element is written -- hence recording the writes rather
  // than applying them as they are found.
  struct ElementChange {
    CoinBigIndex position;
    double value;
  };
  std::vector< ElementChange > elementChanges;

  std::vector< RowBoundChange > rowBoundChanges;
  int firstRow = -1, firstColumn = -1;
  double firstValue = 0.0;

  for (int iRow = 0; iRow < numberRows; iRow++) {
    const CoinBigIndex start = rowStart[iRow];
    const int length = rowLength[iRow];
    // A singleton row carries no slack argument: shrinking its only
    // coefficient would just be a bound change, which is bound propagation's
    // job (and it has already run).
    if (length <= 1)
      continue;

    // Exactly one finite side. A range row has two slacks and neither of them
    // alone justifies the change; a free row has none.
    const bool hasUpper = rowUpper[iRow] < COEFSTR_INFINITY;
    const bool hasLower = rowLower[iRow] > -COEFSTR_INFINITY;
    if (hasUpper == hasLower)
      continue;

    // Normalise to a.x <= rhs: a >= row is negated wholesale.
    const double scale = hasUpper ? 1.0 : -1.0;
    double rhs = hasUpper ? rowUpper[iRow] : -rowLower[iRow];

    // Maximum activity. An infinite contribution makes the row unbounded
    // above, so there is no slack to work with at all.
    double maxActivity = 0.0;
    bool infiniteActivity = false;
    for (CoinBigIndex j = start; j < start + length; j++) {
      const double value = element[j] * scale;
      const int iColumn = column[j];
      if (value > 0.0) {
        if (colUpper[iColumn] >= COEFSTR_INFINITY) {
          infiniteActivity = true;
          break;
        }
        maxActivity += value * colUpper[iColumn];
      } else if (value < 0.0) {
        if (colLower[iColumn] <= -COEFSTR_INFINITY) {
          infiniteActivity = true;
          break;
        }
        maxActivity += value * colLower[iColumn];
      }
    }
    if (infiniteActivity)
      continue;

    // The slack is what every coefficient in this row is capped at. It is
    // invariant under the change itself (the coefficient loses exactly what
    // the rhs loses), so it is computed once and the order of the changes
    // below does not matter.
    double slack = maxActivity - rhs;
    if (slack <= COEFSTR_FEASTOL)
      continue;

    // Snap up to an integral slack when we are within tolerance of one, so a
    // slack that is really 3 but computed as 2.9999995 does not shave a
    // sliver off the feasible set. Rounding is only ever *upwards*: a larger
    // cap is a weaker (safe) tightening, a smaller one can cut off integer
    // points. Non-integral slacks are left alone -- the rule does not need
    // them to be integral.
    const double slackCeil = std::ceil(slack);
    if (std::fabs(slackCeil - slack) <= COEFSTR_FEASTOL)
      slack = slackCeil;

    bool changedRow = false;
    for (CoinBigIndex j = start; j < start + length; j++) {
      const int iColumn = column[j];
      if (!solver->isInteger(iColumn))
        continue;
      // Both bounds must be finite: the rhs adjustment is stated in terms of
      // the bound the column attains at maximum activity.
      if (colLower[iColumn] <= -COEFSTR_INFINITY || colUpper[iColumn] >= COEFSTR_INFINITY)
        continue;

      // The adjustment is stated at an *attainable integer* bound: the
      // argument is that one step in from it the row is already redundant,
      // and a step is only a step if the bound is integral. Bound
      // propagation normally leaves integral bounds on integer columns, but a
      // bound of 0.5 would otherwise give a step of 0.5 and shave a sliver
      // off the feasible set, so round inwards and require a real step to
      // exist -- which also excludes fixed columns, where shrinking the
      // coefficient buys nothing.
      const double lower = std::ceil(colLower[iColumn] - COEFSTR_FEASTOL);
      const double upper = std::floor(colUpper[iColumn] + COEFSTR_FEASTOL);
      if (upper - lower < 1.0 - COEFSTR_FEASTOL)
        continue;

      const double value = element[j] * scale;
      if (std::fabs(value) <= slack + COEFSTR_FEASTOL)
        continue;

      const double newValue = (value > 0.0 ? slack : -slack);
      // At the bound the column takes at maximum activity the two rows agree
      // exactly; one integer step away from it both are redundant.
      rhs -= (value - newValue) * (value > 0.0 ? upper : lower);
      // Recorded, not written: `element` still points at the solver's own row
      // copy. Safe to defer because the slack is invariant under the change
      // and each coefficient is considered exactly once, so nothing later in
      // the scan depends on this write having landed.
      ElementChange elementChange;
      elementChange.position = j;
      elementChange.value = newValue * scale;
      elementChanges.push_back(elementChange);

      if (logLevel >= 3)
        printf("  Coefficient strengthening: row %d col %d: %g -> %g "
               "(slack %g, rhs now %g)\n",
          iRow, iColumn, value, newValue, slack, rhs);

      if (firstRow < 0) {
        firstRow = iRow;
        firstColumn = iColumn;
        firstValue = elementChange.value;
      }
      nCoefficients_++;
      changedRow = true;
    }

    if (changedRow) {
      nRows_++;
      RowBoundChange change;
      change.row = iRow;
      change.lower = hasUpper ? rowLower[iRow] : -rhs;
      change.upper = hasUpper ? rhs : rowUpper[iRow];
      rowBoundChanges.push_back(change);
    }
  }

  if (!nCoefficients_) {
    timeUsed_ = CoinWallclockTime() - t0;
    return false;
  }

  // Only now is the copy worth its cost. Poking the single entries with
  // modifyCoefficient() instead would leave the solver's cached row copy stale
  // -- ClpModel::modifyCoefficient() only clears whatsChanged_ and never calls
  // freeCachedResults() -- so a coefficient written that way reads back with
  // its old value through getMatrixByRow(). replaceMatrix() does free the
  // caches, which is also why Cgl's own tightener batches this way.
  CoinPackedMatrix copy = *rowCopy;
  double *mutableElement = copy.getMutableElements();
  for (size_t i = 0; i < elementChanges.size(); i++)
    mutableElement[elementChanges[i].position] = elementChanges[i].value;

  solver->replaceMatrixOptional(copy);

  // replaceMatrixOptional() is a no-op in the OsiSolverInterface base class,
  // so for a solver that does not implement it the matrix would silently keep
  // its original coefficients. Adjusting the row bounds anyway would then cut
  // off feasible points, so confirm the replacement landed before touching
  // them, and leave the model exactly as it was if it did not.
  bool applied = false;
  {
    const CoinPackedMatrix *check = solver->getMatrixByRow();
    if (check && firstRow >= 0 && firstRow < check->getNumRows()) {
      const CoinBigIndex start = check->getVectorStarts()[firstRow];
      const int length = check->getVectorLengths()[firstRow];
      const int *indices = check->getIndices();
      const double *values = check->getElements();
      for (CoinBigIndex j = start; j < start + length; j++) {
        if (indices[j] == firstColumn) {
          applied = (std::fabs(values[j] - firstValue) < COEFSTR_EPSILON);
          break;
        }
      }
    }
  }

  if (!applied) {
    if (logLevel >= 1)
      printf("  Coefficient strengthening: solver did not accept the modified "
             "matrix — %d changes discarded.\n",
        nCoefficients_);
    nCoefficients_ = 0;
    nRows_ = 0;
    timeUsed_ = CoinWallclockTime() - t0;
    return false;
  }

  for (size_t i = 0; i < rowBoundChanges.size(); i++) {
    const RowBoundChange &change = rowBoundChanges[i];
    if (change.lower > -COEFSTR_INFINITY)
      solver->setRowLower(change.row, change.lower);
    if (change.upper < COEFSTR_INFINITY)
      solver->setRowUpper(change.row, change.upper);
  }

  // Debugger validation: every changed row must still contain the reference
  // solution once evaluated with its NEW coefficients against its NEW
  // bounds. This is a purely algebraic invariant of the rule above (see the
  // class comment) -- unlike a cut, nothing here depends on which subtree of
  // the B&B tree we are in, so it is checked unconditionally whenever a
  // debugger/-debugCuts solution is available, not just on the optimal path.
  //
  // Fall back to the debugSolution global (see the comment above the
  // extern declarations) when no OsiRowCutDebugger has been attached to
  // `solver` yet -- this pass runs ahead of the root LP solve, so a
  // "-debugCuts" driven debugger is not attached this early, while
  // Cbc_activateRowCutDebugger()/mip-debug-cuts attach it before Cbc_solve()
  // is even called, so a debugger already exists there.
  //
  // Deliberately use getRowCutDebuggerAlways() and NOT getRowCutDebugger():
  // the latter gates on OsiRowCutDebugger::onOptimalPath(), which only
  // checks INTEGER columns with a loose 1e-3 tolerance and has been observed
  // to be unreliable (can return false even though the current bounds are
  // still fully consistent with the reference solution) -- see the matching
  // fixes/comments in CbcModel.cpp (reducedCostFix(), fixFromGlobalCuts())
  // and CbcBoundPropagation::run(). Using the unreliable check here would
  // silently skip this validation whenever it spuriously returns false,
  // which is exactly the failure mode this diagnostic exists to catch.
  const OsiRowCutDebugger *debugger = solver->getRowCutDebuggerAlways();
  const double *optSol = debugger
    ? debugger->optimalSolution()
    : (debugSolution && debugNumberColumns == solver->getNumCols()
         ? debugSolution
         : nullptr);
  if (optSol) {
    const double *checkElement = copy.getElements();
    const int *checkColumn = copy.getIndices();
    const CoinBigIndex *checkRowStart = copy.getVectorStarts();
    const int *checkRowLength = copy.getVectorLengths();
    const double tol = 1.0e-6;
    for (size_t i = 0; i < rowBoundChanges.size(); i++) {
      const RowBoundChange &change = rowBoundChanges[i];
      const int iRow = change.row;
      double activity = 0.0;
      const CoinBigIndex start = checkRowStart[iRow];
      const int length = checkRowLength[iRow];
      for (CoinBigIndex j = start; j < start + length; j++)
        activity += checkElement[j] * optSol[checkColumn[j]];
      const bool violatesUpper = change.upper < COEFSTR_INFINITY && activity > change.upper + tol;
      const bool violatesLower = change.lower > -COEFSTR_INFINITY && activity < change.lower - tol;
      if (violatesUpper || violatesLower) {
        printf("coefStrengthening BAD ROW: row %d activity=%.12g at reference "
               "solution but new bounds are [%.12g,%.12g] (violates %s)\n",
          iRow, activity, change.lower, change.upper,
          violatesUpper ? "upper" : "lower");
        fflush(stdout);
      }
    }
  }

  timeUsed_ = CoinWallclockTime() - t0;
  return true;
}
