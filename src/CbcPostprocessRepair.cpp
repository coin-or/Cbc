// CbcPostprocessRepair.cpp
//
// Implementation of the three-phase postprocess repair pass.
// See CbcPostprocessRepair.hpp for a high-level description.

#include "CbcPostprocessRepair.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "CbcModel.hpp"
#include "CglPreProcess.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"

void CbcRepairPostprocessSolution(
  OsiSolverInterface *saveSolver,
  OsiSolverInterface *originalSolver,
  CbcModel *babModel,
  CglPreProcess &process)
{
  int nCols = originalSolver->getNumCols();
  int nColsSave = saveSolver->getNumCols();
  // Only run the repair if saveSolver has at least as many columns as
  // originalSolver (= original problem space).  nColsSave < nCols would
  // mean postProcess changed the model unexpectedly.
  if (nColsSave < nCols)
    return;

  // Start from the postProcess back-substituted solution.
  // postProcess correctly handles ALL preprocessing transformations including
  // variable sign-flips (complementing), so we trust it for all variables.
  const double *proposedSol = saveSolver->getColSolution();
  std::vector< double > repSol(proposedSol, proposedSol + nColsSave);

  const CoinPackedMatrix *rowMtx = originalSolver->getMatrixByRow();
  const CoinPackedMatrix *colMtx = originalSolver->getMatrixByCol();
  const double *rowLb = originalSolver->getRowLower();
  const double *rowUb = originalSolver->getRowUpper();
  const double *origColLb = originalSolver->getColLower();
  const double *origColUb = originalSolver->getColUpper();
  int nRows = originalSolver->getNumRows();
  std::vector< double > lhs(nRows, 0.0);
  for (int r = 0; r < nRows; r++) {
    const int s = rowMtx->getVectorStarts()[r];
    const int rl = rowMtx->getVectorLengths()[r];
    for (int k = 0; k < rl; k++)
      lhs[r] += rowMtx->getElements()[s + k] * repSol[rowMtx->getIndices()[s + k]];
  }

  // locked[col]: variable committed; do not change again.
  std::vector< bool > locked(nCols, false);
  std::vector< bool > changed(nCols, false);

  // Pre-compute pureInt[r]: true if every variable in row r is integer.
  // Only pure-integer rows are repaired; rows with continuous variables
  // are handled by the subsequent LP re-solve.
  std::vector< bool > pureInt(nRows, true);
  for (int r = 0; r < nRows; r++) {
    const int s = rowMtx->getVectorStarts()[r];
    const int rl = rowMtx->getVectorLengths()[r];
    for (int k = 0; k < rl; k++) {
      if (!originalSolver->isInteger(rowMtx->getIndices()[s + k])) {
        pureInt[r] = false;
        break;
      }
    }
  }

  // Pre-compute forced values from singleton rows in the original model.
  // A singleton row (exactly 1 variable, coefficient != 0) with a tight RHS
  // uniquely determines that variable's value from the constraint alone.
  // These forced values take priority over oracle values, because OsiPresolve
  // may have COMPLEMENTED variables (transforming x → 1-x or x → ub-x) so
  // that babBest[pre] is the VALUE OF THE COMPLEMENT, not the original variable.
  // Singleton rows always encode the original-model constraint faithfully.
  std::vector< double > singletonForced(nCols, std::numeric_limits< double >::quiet_NaN());
  for (int r = 0; r < nRows; r++) {
    const int rs = rowMtx->getVectorStarts()[r];
    const int rl = rowMtx->getVectorLengths()[r];
    if (rl != 1)
      continue;
    int col = rowMtx->getIndices()[rs];
    if (!originalSolver->isInteger(col))
      continue;
    double coeff = rowMtx->getElements()[rs];
    if (fabs(coeff) < 1e-12)
      continue;
    double rLb = rowLb[r], rUb = rowUb[r];
    // Range for col value implied by the row constraint
    double colRangeMin = (coeff > 0) ? rLb / coeff : rUb / coeff;
    double colRangeMax = (coeff > 0) ? rUb / coeff : rLb / coeff;
    // Intersect with variable bounds
    colRangeMin = std::max(colRangeMin, origColLb[col]);
    colRangeMax = std::min(colRangeMax, origColUb[col]);
    double lo = ceil(colRangeMin - 1e-8);
    double hi = floor(colRangeMax + 1e-8);
    if (lo > hi + 0.5)
      continue; // infeasible singleton — skip
    if (hi - lo < 0.5) { // unique integer forced
      singletonForced[col] = lo;
    }
  }

  // Build an oracle: for each original-model integer column, the correct
  // integer value as determined by the B&B best solution.
  // CglPreProcess::originalColumns() maps preprocessed-model column index →
  // original-model column index, so we invert it here.
  // babModel->bestSolution() is in the presolved (OsiPresolve) variable space.
  // We apply oracle to ALL presolved-mapped variables — not just fractional ones —
  // because the LP re-solves in postProcess can also leave non-fractional integer
  // variables at LP-optimal values that differ from the B&B optimal.  Applying
  // all 1100 oracle anchors before propagation gives the correct starting point
  // for reconstructing the full original-space B&B solution.
  // NOTE: if a singleton row forces a different value than the oracle, the singleton
  // wins (oracle can be incorrect for OsiPresolve-complemented variables).
  std::vector< double > oracle(nCols, std::numeric_limits< double >::quiet_NaN());
  if (babModel && babModel->bestSolution()) {
    const double *babBest = babModel->bestSolution();
    int nBabCols = babModel->getNumCols();
    const int *origColsMap = process.originalColumns(); // preprocessed → original
    for (int pre = 0; pre < nBabCols; pre++) {
      int orig = origColsMap ? origColsMap[pre] : pre;
      if (orig >= 0 && orig < nCols && originalSolver->isInteger(orig))
        oracle[orig] = floor(babBest[pre] + 0.5);
    }
  }

  // anchored[col]: variable has been set to a trusted integer value (either
  // directly from the B&B oracle, fixed by variable bounds lb==ub, or forced
  // by a singleton row constraint).
  // Anchored variables are used as the base for constraint propagation to
  // reconstruct the values of variables that were eliminated by OsiPresolve.
  std::vector< bool > anchored(nCols, false);

  // Phase 1: set integer variables to their B&B optimal or constraint-forced
  // values. Priority order:
  //   1. Bounds-fixed (lb == ub): variable is fixed regardless.
  //   2. Singleton-forced: a single-variable row uniquely determines the value.
  //      This OVERRIDES oracle because OsiPresolve may have complemented the
  //      variable (oracle gives the complement value, not the original).
  //   3. Oracle: B&B solution value from presolved space.
  // Variables with none of the above are left for Phase 1b propagation.
  {
    for (int col = 0; col < nCols; col++) {
      if (!originalSolver->isInteger(col))
        continue;
      double lb = origColLb[col];
      double ub = origColUb[col];
      if (ub - lb < 0.5) {
        // Bounds-fixed variable: anchor at midpoint.
        double fixedVal = floor(lb + 0.5);
        double delta = fixedVal - repSol[col];
        if (fabs(delta) > 1e-10) {
          const int cs = colMtx->getVectorStarts()[col];
          const int cl = colMtx->getVectorLengths()[col];
          for (int k = 0; k < cl; k++)
            lhs[colMtx->getIndices()[cs + k]] += delta * colMtx->getElements()[cs + k];
          repSol[col] = fixedVal;
          changed[col] = true;
        }
        anchored[col] = true;
        continue;
      }
      // Determine the value to anchor: singleton-forced takes priority over oracle.
      double anchorVal = std::numeric_limits< double >::quiet_NaN();
      bool lockIt = false;
      if (!std::isnan(singletonForced[col])) {
        anchorVal = singletonForced[col]; // constraint-forced value
        lockIt = true; // singleton-pinned: Phase 2 must not move this variable
      } else if (!std::isnan(oracle[col])) {
        anchorVal = floor(oracle[col] + 0.5);
        anchorVal = std::max(lb, std::min(ub, anchorVal));
      }
      if (std::isnan(anchorVal))
        continue; // no oracle and no singleton force — handle in Phase 1b
      double delta = anchorVal - repSol[col];
      if (fabs(delta) > 1e-10) {
        const int cs = colMtx->getVectorStarts()[col];
        const int cl = colMtx->getVectorLengths()[col];
        for (int k = 0; k < cl; k++)
          lhs[colMtx->getIndices()[cs + k]] += delta * colMtx->getElements()[cs + k];
        repSol[col] = anchorVal;
        changed[col] = true;
      }
      anchored[col] = true;
      if (lockIt)
        locked[col] = true;
    }
  }

  // Phase 1.5: Conflict propagation from locked (singleton-constrained) variables.
  // Oracle values from the presolved space may be wrong for OsiPresolve-complemented
  // variables.  When a locked variable forces a row constraint to be infeasible
  // (because the oracle set another variable to the wrong value), we adjust the
  // non-locked variable to restore feasibility.
  // This is the "bound propagation" step the user described: if fixing locked vars
  // + oracle vars creates an infeasibility, it definitely means some oracle value
  // is wrong, so we correct it here.
  // Iterates until no more corrections are possible.
  {
    bool progress = true;
    while (progress) {
      progress = false;
      for (int r = 0; r < nRows; r++) {
        if (!pureInt[r])
          continue;
        double viol_lo = rowLb[r] - lhs[r]; // positive: lhs too low
        double viol_hi = lhs[r] - rowUb[r]; // positive: lhs too high
        if (viol_lo <= 1e-6 && viol_hi <= 1e-6)
          continue; // row satisfied

        const int rs = rowMtx->getVectorStarts()[r];
        const int rl = rowMtx->getVectorLengths()[r];

        // Only adjust if row contains at least one locked variable
        // (i.e., the violation is caused by a locked var constraining others).
        bool hasLocked = false;
        for (int k = 0; k < rl; k++)
          if (locked[rowMtx->getIndices()[rs + k]]) {
            hasLocked = true;
            break;
          }
        if (!hasLocked)
          continue;

        // Correct non-locked integer variables that contribute to the violation.
        for (int k = 0; k < rl; k++) {
          int c = rowMtx->getIndices()[rs + k];
          if (locked[c])
            continue;
          if (!originalSolver->isInteger(c))
            continue;
          double coeff = rowMtx->getElements()[rs + k];
          if (fabs(coeff) < 1e-12)
            continue;

          double cur = repSol[c];
          double cLb = origColLb[c];
          double cUb = origColUb[c];
          double newVal = cur;

          if (viol_hi > 1e-6 && coeff > 0 && cur > cLb + 1e-8) {
            // lhs too high: decrease this variable toward its lower bound.
            double decrease = std::min(viol_hi / coeff, cur - cLb);
            newVal = std::max(cLb, floor(cur - decrease + 1e-8));
          } else if (viol_lo > 1e-6 && coeff > 0 && cur < cUb - 1e-8) {
            // lhs too low: increase this variable toward its upper bound.
            double increase = std::min(viol_lo / coeff, cUb - cur);
            newVal = std::min(cUb, ceil(cur + increase - 1e-8));
          } else if (viol_hi > 1e-6 && coeff < 0 && cur < cUb - 1e-8) {
            // lhs too high: increase negative-coeff var to decrease lhs.
            double increase = std::min(viol_hi / (-coeff), cUb - cur);
            newVal = std::min(cUb, ceil(cur + increase - 1e-8));
          } else if (viol_lo > 1e-6 && coeff < 0 && cur > cLb + 1e-8) {
            // lhs too low: decrease negative-coeff var to increase lhs.
            double decrease = std::min(viol_lo / (-coeff), cur - cLb);
            newVal = std::max(cLb, floor(cur - decrease + 1e-8));
          } else {
            continue;
          }

          double delta = newVal - cur;
          if (fabs(delta) < 1e-10)
            continue;
          const int cs = colMtx->getVectorStarts()[c];
          const int cl = colMtx->getVectorLengths()[c];
          for (int kk = 0; kk < cl; kk++)
            lhs[colMtx->getIndices()[cs + kk]] += delta * colMtx->getElements()[cs + kk];
          repSol[c] = newVal;
          changed[c] = true;
          anchored[c] = true; // re-anchor at corrected value
          // Recompute violations for this row after change
          viol_lo = rowLb[r] - lhs[r];
          viol_hi = lhs[r] - rowUb[r];
          progress = true;
        }
      }
    }
  }

  // Phase 1b: constraint propagation for non-anchored integer variables.
  // Variables eliminated by OsiPresolve have no oracle.  Many appear in
  // pure-integer rows where all OTHER variables are now anchored.  Scan
  // iteratively for such "singleton" rows and fix the remaining variable.
  // This reconstructs the back-substitution that OsiPresolve would perform,
  // without needing access to the presolve chain.
  // Handles set-partitioning/covering patterns: sum(x_i)=1, where oracle
  // anchored k-1 variables and the remaining one is uniquely determined.
  {
    bool progress = true;
    while (progress) {
      progress = false;
      for (int r = 0; r < nRows; r++) {
        if (!pureInt[r])
          continue;
        const int rs = rowMtx->getVectorStarts()[r];
        const int rl = rowMtx->getVectorLengths()[r];
        // Find the single non-anchored variable in this row.
        int freeCol = -1;
        int freeCount = 0;
        double residual = 0.0;
        for (int k = 0; k < rl; k++) {
          int col = rowMtx->getIndices()[rs + k];
          if (!anchored[col]) {
            freeCount++;
            freeCol = col;
          } else {
            residual += rowMtx->getElements()[rs + k] * repSol[col];
          }
        }
        if (freeCount != 1)
          continue;

        double coeff = 0.0;
        for (int k = 0; k < rl; k++)
          if (rowMtx->getIndices()[rs + k] == freeCol)
            coeff = rowMtx->getElements()[rs + k];
        if (fabs(coeff) < 1e-12)
          continue;

        double lb2 = origColLb[freeCol];
        double ub2 = origColUb[freeCol];
        double requiredLo = (rowLb[r] - residual) / coeff;
        double requiredHi = (rowUb[r] - residual) / coeff;
        if (coeff < 0)
          std::swap(requiredLo, requiredHi);
        requiredLo = std::max(lb2, ceil(requiredLo - 1e-8));
        requiredHi = std::min(ub2, floor(requiredHi + 1e-8));
        if (requiredLo > requiredHi + 0.5)
          continue;

        double nearest = floor(repSol[freeCol] + 0.5);
        double chosen = std::max(requiredLo, std::min(requiredHi, nearest));
        double delta = chosen - repSol[freeCol];
        if (fabs(delta) > 1e-10) {
          const int cs = colMtx->getVectorStarts()[freeCol];
          const int cl = colMtx->getVectorLengths()[freeCol];
          for (int k = 0; k < cl; k++)
            lhs[colMtx->getIndices()[cs + k]] += delta * colMtx->getElements()[cs + k];
          repSol[freeCol] = chosen;
          changed[freeCol] = true;
        }
        anchored[freeCol] = true;
        progress = true;
      }
    }
  }

  // Phase 1c: nearest-integer fallback for any remaining non-anchored integer
  // variables that constraint propagation could not determine.  Phase 2/3 will
  // repair any constraint violations introduced by this rounding.
  {
    for (int col = 0; col < nCols; col++) {
      if (!originalSolver->isInteger(col) || anchored[col])
        continue;
      double v = repSol[col];
      double fl = floor(v);
      double ce = fl + 1.0;
      double lb = origColLb[col];
      double ub = origColUb[col];
      if (ub - lb < 0.5)
        continue;
      double rounded = (v - fl >= 0.5) ? ce : fl;
      rounded = std::max(lb, std::min(ub, rounded));
      double delta = rounded - v;
      if (fabs(delta) < 1e-10)
        continue;
      const int cs = colMtx->getVectorStarts()[col];
      const int cl = colMtx->getVectorLengths()[col];
      for (int k = 0; k < cl; k++)
        lhs[colMtx->getIndices()[cs + k]] += delta * colMtx->getElements()[cs + k];
      repSol[col] = rounded;
      changed[col] = true;
      anchored[col] = true;
    }
  }

  // Phase 2: variable-first greedy repair (up to nCols passes).
  // Each pass picks the integer variable whose change reduces total
  // row-violation the most across pure-integer rows, commits it,
  // and repeats.  Handles structural violations that remain after
  // Phase 1 rounding (e.g., binary variables still at the wrong
  // endpoint, or coupled integer variables that rounding alone
  // could not resolve).
  for (int pass = 0; pass < 2 * nCols; pass++) {
    // Check if any pure-integer rows are still violated.
    bool anyViolated = false;
    for (int r = 0; r < nRows; r++) {
      if (!pureInt[r])
        continue;
      if (std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]) > 1.0e-6) {
        anyViolated = true;
        break;
      }
    }
    if (!anyViolated)
      break;

    int bestCol = -1;
    double bestGain = 1.0e-8; // must strictly improve
    double bestNewVal = 0.0;

    for (int col = 0; col < nCols; col++) {
      if (locked[col])
        continue;
      if (!originalSolver->isInteger(col))
        continue;
      double curVal = repSol[col];
      double colLb = origColLb[col];
      double colUb = origColUb[col];
      double range = colUb - colLb;
      if (range < 0.5)
        continue; // fixed variable, skip

      const int cs = colMtx->getVectorStarts()[col];
      const int cl = colMtx->getVectorLengths()[col];
      const int *cIdx = colMtx->getIndices() + cs;
      const double *cElem = colMtx->getElements() + cs;

      // Build candidate new values for this variable.
      // For binaries: the opposite endpoint.
      // For general integers: projection-based ideal values (one
      // per violated row containing this column) and ±1/±2 steps.
      // Using a small fixed-size buffer avoids heap allocation in
      // the inner loop.
      double candBuf[16];
      int nCands = 0;

      if (range <= 1.5) {
        // Binary: try the other 0/1 endpoint.
        candBuf[nCands++] = (curVal < colLb + 0.5) ? (colLb + 1.0) : colLb;
      } else {
        // General integer: collect projection candidates from each
        // violated pure-integer row that contains this column.
        for (int k = 0; k < cl && nCands < 12; k++) {
          int r = cIdx[k];
          if (!pureInt[r])
            continue;
          double coeff = cElem[k];
          if (fabs(coeff) < 1.0e-10)
            continue;
          double viol = 0.0;
          if (lhs[r] < rowLb[r] - 1.0e-6)
            viol = rowLb[r] - lhs[r]; // positive: lhs too small
          else if (lhs[r] > rowUb[r] + 1.0e-6)
            viol = rowUb[r] - lhs[r]; // negative: lhs too large
          else
            continue; // row not violated
          // Ideal delta: viol / coeff moves this row exactly to bound.
          double idealVal = curVal + viol / coeff;
          double fl = floor(idealVal + 1.0e-8);
          double ce = ceil(idealVal - 1.0e-8);
          for (double cand : { fl, ce }) {
            cand = std::max(colLb, std::min(colUb, cand));
            if (fabs(cand - curVal) > 1.0e-8 && nCands < 14)
              candBuf[nCands++] = cand;
          }
        }
        // Also try ±1 and ±2 for small back-substitution errors.
        for (double delta : { -1.0, +1.0, -2.0, +2.0 }) {
          double cand = std::max(colLb, std::min(colUb, curVal + delta));
          if (fabs(cand - curVal) > 1.0e-8 && nCands < 16)
            candBuf[nCands++] = cand;
        }
      }

      // Evaluate each candidate and track the best gain.
      for (int ci = 0; ci < nCands; ci++) {
        double newVal = candBuf[ci];
        double delta = newVal - curVal;
        if (fabs(delta) < 1.0e-8)
          continue;
        double gain = 0.0;
        for (int k = 0; k < cl; k++) {
          int r = cIdx[k];
          if (!pureInt[r])
            continue;
          double curViol = std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]);
          if (curViol < 0.0)
            curViol = 0.0;
          double newLhs = lhs[r] + delta * cElem[k];
          double newViol = std::max(rowLb[r] - newLhs, newLhs - rowUb[r]);
          if (newViol < 0.0)
            newViol = 0.0;
          gain += curViol - newViol;
        }
        if (gain > bestGain) {
          bestGain = gain;
          bestCol = col;
          bestNewVal = newVal;
        }
      }
    }

    if (bestCol < 0)
      break; // no improving move found

    // Commit the best move: update lhs for affected rows.
    double delta = bestNewVal - repSol[bestCol];
    const int cs = colMtx->getVectorStarts()[bestCol];
    const int cl = colMtx->getVectorLengths()[bestCol];
    const int *cIdx = colMtx->getIndices() + cs;
    const double *cElem = colMtx->getElements() + cs;
    for (int k = 0; k < cl; k++)
      lhs[cIdx[k]] += delta * cElem[k];
    repSol[bestCol] = bestNewVal;
    locked[bestCol] = true;
    changed[bestCol] = true;
  }

  // Phase 3: Tabu search for remaining violations.
  //
  // After Phase 1+2, remaining violations have gain=0 for every single
  // variable move (fix one row, break another).  Tabu search escapes this
  // local minimum by allowing zero/negative gain moves while preventing
  // cycling via a tabu list that blocks reversal for tabuTenure iterations.
  //
  // Only variables appearing in currently-violated rows are considered as
  // candidates, keeping per-iteration cost O(violRows × rowLen × colDensity).
  {
    int nViolP3start = 0;
    for (int r = 0; r < nRows; r++)
      if (pureInt[r]
        && std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]) > 1e-6)
        nViolP3start++;

    if (nViolP3start > 0) {
      // Save current state so we can restore the best-seen solution.
      std::vector< double > repSolP3 = repSol;
      std::vector< double > bestRepSol = repSol;
      std::vector< double > bestLhs = lhs;
      int bestViol = nViolP3start;

      std::vector< int > tabuEnd(nCols, -1);
      const int tabuTenure = 7;
      const int maxIter = 10000;

      int curViol = nViolP3start;

      for (int iter = 0; iter < maxIter; iter++) {
        if (curViol == 0)
          break;

        // Scan violated rows to collect candidate columns.
        int bestCol = -1, bestDir = 0;
        double bestGain = -1e30;
        std::vector< bool > seen(nCols, false);

        for (int r = 0; r < nRows; r++) {
          if (!pureInt[r])
            continue;
          if (std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]) <= 1e-6)
            continue;
          const int rs = rowMtx->getVectorStarts()[r];
          const int rl = rowMtx->getVectorLengths()[r];
          for (int ki = 0; ki < rl; ki++) {
            int c = rowMtx->getIndices()[rs + ki];
            if (seen[c])
              continue;
            seen[c] = true;
            if (!originalSolver->isInteger(c))
              continue;
            if (locked[c])
              continue;

            for (int dir : { -1, 1 }) {
              double nv = repSol[c] + dir;
              if (nv < origColLb[c] - 1e-10 || nv > origColUb[c] + 1e-10)
                continue;

              // Compute gain (violation reduction).
              double gain = 0.0;
              const int cs = colMtx->getVectorStarts()[c];
              const int cl = colMtx->getVectorLengths()[c];
              for (int kk = 0; kk < cl; kk++) {
                int rr = colMtx->getIndices()[cs + kk];
                if (!pureInt[rr])
                  continue;
                double v0 = std::max(0.0, std::max(rowLb[rr] - lhs[rr], lhs[rr] - rowUb[rr]));
                double newLhs = lhs[rr] + dir * colMtx->getElements()[cs + kk];
                double v1 = std::max(0.0, std::max(rowLb[rr] - newLhs, newLhs - rowUb[rr]));
                gain += v0 - v1;
              }

              bool isTabu = (tabuEnd[c] > iter);
              // Aspiration: accept tabu move if it beats best solution.
              bool aspiration = (curViol - gain) < bestViol;
              if (gain > bestGain && (!isTabu || aspiration)) {
                bestGain = gain;
                bestCol = c;
                bestDir = dir;
              }
            }
          }
        }

        if (bestCol == -1)
          break;

        // Apply best move.
        const int cs = colMtx->getVectorStarts()[bestCol];
        const int cl = colMtx->getVectorLengths()[bestCol];
        for (int kk = 0; kk < cl; kk++)
          lhs[colMtx->getIndices()[cs + kk]] += bestDir * colMtx->getElements()[cs + kk];
        repSol[bestCol] += bestDir;
        tabuEnd[bestCol] = iter + tabuTenure;

        // Recount violations.
        curViol = 0;
        for (int r = 0; r < nRows; r++)
          if (pureInt[r] && std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]) > 1e-6)
            curViol++;

        if (curViol < bestViol) {
          bestViol = curViol;
          bestRepSol = repSol;
          bestLhs = lhs;
        }
      }

      // Restore best solution found during tabu search.
      repSol = bestRepSol;
      lhs = bestLhs;

      // Phase 3b: depth-limited DFS repair.
      // Allows violations to temporarily increase by up to 2 beyond bestViol
      // to escape local minima that tabu search could not resolve.
      // Uses incremental violation tracking with a violated-row set to keep
      // each DFS node cost proportional to column density (not nRows).
      if (bestViol > 0 && bestViol <= 6) {
        std::vector< bool > isViol3b(nRows, false);
        std::vector< int > violRows3b;
        for (int r = 0; r < nRows; r++) {
          if (pureInt[r]
            && std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]) > 1e-6) {
            isViol3b[r] = true;
            violRows3b.push_back(r);
          }
        }

        // Apply move (c,d) with incremental violated-row update.
        // Returns delta (change in violation count).
        auto apply3b = [&](int c, int d) -> int {
          repSol[c] += d;
          int delta = 0;
          const int cs2 = colMtx->getVectorStarts()[c];
          const int cl2 = colMtx->getVectorLengths()[c];
          for (int kk2 = 0; kk2 < cl2; kk2++) {
            int rr = colMtx->getIndices()[cs2 + kk2];
            if (!pureInt[rr])
              continue;
            bool wasV = isViol3b[rr];
            lhs[rr] += d * colMtx->getElements()[cs2 + kk2];
            bool isV = std::max(rowLb[rr] - lhs[rr], lhs[rr] - rowUb[rr]) > 1e-6;
            if (wasV && !isV) {
              isViol3b[rr] = false;
              violRows3b.erase(std::find(violRows3b.begin(), violRows3b.end(), rr));
              delta--;
            } else if (!wasV && isV) {
              isViol3b[rr] = true;
              violRows3b.push_back(rr);
              delta++;
            }
          }
          return delta;
        };
        // Revert: just apply the reverse.
        auto revert3b = [&](int c, int d) { apply3b(c, -d); };

        const int dfsMaxDepth = 8;
        const int dfsMaxViol = bestViol + 2;
        const int dfsNodeLimit = 50000;
        int dfsNodes = 0;
        bool solved3b = false;

        std::function< bool(int, int, int, int) > dfs3b = [&](int depth, int curV, int prevC, int prevD) -> bool {
          if (curV == 0)
            return true;
          if (depth >= dfsMaxDepth || dfsNodes >= dfsNodeLimit)
            return false;
          // Collect candidate (col,dir) pairs from violated rows.
          std::vector< int > cCols;
          for (int r : violRows3b) {
            const int rs3 = rowMtx->getVectorStarts()[r];
            const int rl3 = rowMtx->getVectorLengths()[r];
            for (int ki3 = 0; ki3 < rl3; ki3++) {
              int cc = rowMtx->getIndices()[rs3 + ki3];
              if (originalSolver->isInteger(cc) && !locked[cc])
                cCols.push_back(cc);
            }
          }
          std::sort(cCols.begin(), cCols.end());
          cCols.erase(std::unique(cCols.begin(), cCols.end()), cCols.end());
          for (int c : cCols) {
            for (int d : { -1, 1 }) {
              if (c == prevC && d == -prevD)
                continue;
              double nval = repSol[c] + d;
              if (nval < origColLb[c] - 1e-10 || nval > origColUb[c] + 1e-10)
                continue;
              ++dfsNodes;
              int delta = apply3b(c, d);
              int newV = curV + delta;
              if (newV <= dfsMaxViol && dfs3b(depth + 1, newV, c, d)) {
                return true;
              }
              revert3b(c, d);
            }
          }
          return false;
        };

        solved3b = dfs3b(0, bestViol, -1, 0);

        if (solved3b) {
          bestViol = 0;
          bestRepSol = repSol;
          bestLhs = lhs;
        } else {
          repSol = bestRepSol;
          lhs = bestLhs;
        }
      }

      // Update changed[] for any variable that moved from Phase-2 state.
      for (int c = 0; c < nCols; c++) {
        if (fabs(repSol[c] - repSolP3[c]) > 0.5)
          changed[c] = true;
      }
    }
  }

  // Inject the final repSol back into saveSolver and fix integer variable
  // bounds so that subsequent LP solves preserve the repaired values.
  for (int col = 0; col < nCols; col++) {
    if (saveSolver->isInteger(col)) {
      double val = floor(repSol[col] + 0.5);
      saveSolver->setColLower(col, val);
      saveSolver->setColUpper(col, val);
    } else if (changed[col]) {
      // For non-integer changed variables (shouldn't happen in practice
      // since we only flip binaries, but be safe).
      saveSolver->setColLower(col, repSol[col]);
      saveSolver->setColUpper(col, repSol[col]);
    }
  }
  saveSolver->setColSolution(repSol.data());
}
