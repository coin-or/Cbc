// CbcPostprocessRepair.cpp
//
// Implementation of the three-phase postprocess repair pass.
// See CbcPostprocessRepair.hpp for a high-level description.

#include "CbcPostprocessRepair.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
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

  // Build repSol using direct column mapping from the preprocessed
  // B&B solution wherever possible, overriding the back-substituted
  // saveSolver values for retained variables.
  //
  // Background: process.postProcess() back-substitutes the preprocessed
  // solution into the original variable space.  For variables that were
  // ELIMINATED during preprocessing, the back-substitution uses the
  // internally stored (and possibly sign-flipped) constraint equations,
  // which gives wrong 0/1 values when a G-constraint was scaled to an
  // L-constraint.  For RETAINED variables (those that survive into the
  // preprocessed model), the correct values are available directly from
  // babModel->bestSolution() via process.originalColumns().
  //
  // Strategy:
  //  1. Start from saveSolver->getColSolution() (correct for retained).
  //  2. Override retained-variable values with the direct B&B mapping.
  //  3. Run a greedy repair pass only for pure-binary rows.
  const double *proposedSol = saveSolver->getColSolution();
  std::vector< double > repSol(proposedSol, proposedSol + nColsSave);
  // Keep a copy of the postProcess back-substituted values for diagnostics.
  std::vector< double > postProcessSol(repSol);

  // Mark which original-space columns are retained in the preprocessed model.
  std::vector< bool > isRetainedEarly(nCols, false);
  if (babModel->bestSolution()) {
    int nPrep = babModel->solver()->getNumCols();
    const double *prepBest = babModel->bestSolution();
    const int *origCols = process.originalColumns();
    for (int i = 0; i < nPrep; i++) {
      int j = origCols[i]; // original-space column
      if (j < 0 || j >= nColsSave)
        continue;
      isRetainedEarly[j] = true;
      double prepVal = prepBest[i];
      if (babModel->solver()->isInteger(i))
        prepVal = floor(prepVal + 0.5);
      repSol[j] = prepVal;
    }
  }

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
  // With retained variables correctly mapped from the B&B solution, we
  // can trust any row where all variables are integer.  Rows with
  // continuous variables are still excluded (LP will handle those).
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

  // Diagnostic: count violations pre/post retained-variable override.
  // This helps determine if the override is introducing violations.
  {
    // Count pure-int violations in postProcessSol (before override).
    int violPre = 0, violPost = 0;
    for (int r = 0; r < nRows; r++) {
      if (!pureInt[r])
        continue;
      double lhsPre = 0.0, lhsPost = 0.0;
      const int s = rowMtx->getVectorStarts()[r];
      const int rl = rowMtx->getVectorLengths()[r];
      for (int k = 0; k < rl; k++) {
        int c = rowMtx->getIndices()[s + k];
        double e = rowMtx->getElements()[s + k];
        if (c < (int)postProcessSol.size())
          lhsPre += e * postProcessSol[c];
        lhsPost += e * repSol[c];
      }
      if (std::max(rowLb[r] - lhsPre, lhsPre - rowUb[r]) > 1e-6)
        violPre++;
      if (std::max(rowLb[r] - lhsPost, lhsPost - rowUb[r]) > 1e-6)
        violPost++;
    }
    if (violPre > 0 || violPost > 0)
      fprintf(stderr, "[CbcRepairPostprocess] violations pre-override=%d post-override=%d\n",
        violPre, violPost);
  }
  //
  // OsiPresolve back-substitution can leave eliminated integer variables
  // with fractional values (cascading floating-point errors through the
  // multi-pass postsolve chain).  Rounding each such variable in the
  // direction that most reduces the total violation across all rows it
  // appears in fixes the vast majority of violations without any search.
  // Retained variables (from origCols) already have exact B&B values
  // and are skipped.
  {
    // Mark retained variables.
    std::vector< bool > isRetained(nCols, false);
    if (babModel->bestSolution()) {
      int nPrep = babModel->solver()->getNumCols();
      const int *origCols = process.originalColumns();
      for (int i = 0; i < nPrep; i++) {
        int j = origCols[i];
        if (j >= 0 && j < nCols)
          isRetained[j] = true;
      }
    }

    for (int col = 0; col < nCols; col++) {
      if (isRetained[col])
        continue;
      if (!originalSolver->isInteger(col))
        continue;
      double fl = floor(repSol[col]);
      double ce = fl + 1.0;
      // Skip if already integer.
      if (fabs(repSol[col] - fl) < 1e-10 || fabs(repSol[col] - ce) < 1e-10)
        continue;

      double lb = origColLb[col];
      double ub = origColUb[col];
      if (ub - lb < 0.5)
        continue; // fixed variable

      const int cs = colMtx->getVectorStarts()[col];
      const int cl = colMtx->getVectorLengths()[col];
      const int *cIdx = colMtx->getIndices() + cs;
      const double *cElem = colMtx->getElements() + cs;
      double deltaDown = fl - repSol[col]; // < 0
      double deltaUp = ce - repSol[col]; // > 0

      // Compute net violation change for rounding down vs up.
      double gainDown = 0.0, gainUp = 0.0;
      for (int k = 0; k < cl; k++) {
        int r = cIdx[k];
        if (!pureInt[r])
          continue;
        double curViol = std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]);
        if (curViol <= 0.0)
          continue;
        double adr = deltaDown * cElem[k];
        double aur = deltaUp * cElem[k];
        double newLhsD = lhs[r] + adr;
        double newLhsU = lhs[r] + aur;
        gainDown += curViol - std::max(0.0, std::max(rowLb[r] - newLhsD, newLhsD - rowUb[r]));
        gainUp += curViol - std::max(0.0, std::max(rowLb[r] - newLhsU, newLhsU - rowUb[r]));
      }

      double bestRound;
      if (gainDown >= gainUp && fl >= lb)
        bestRound = fl;
      else if (ce <= ub)
        bestRound = ce;
      else if (fl >= lb)
        bestRound = fl;
      else
        continue; // out of bounds both directions — skip

      // Clamp to bounds.
      bestRound = std::max(lb, std::min(ub, bestRound));
      double delta = bestRound - repSol[col];
      if (fabs(delta) < 1e-10)
        continue;

      // Commit the rounding.
      for (int k = 0; k < cl; k++)
        lhs[cIdx[k]] += delta * cElem[k];
      repSol[col] = bestRound;
      changed[col] = true;
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
        bool solved3b = false;

        std::function< bool(int, int, int, int) > dfs3b = [&](int depth, int curV, int prevC, int prevD) -> bool {
          if (curV == 0)
            return true;
          if (depth >= dfsMaxDepth)
            return false;
          // Collect candidate (col,dir) pairs from violated rows.
          std::vector< int > cCols;
          for (int r : violRows3b) {
            const int rs3 = rowMtx->getVectorStarts()[r];
            const int rl3 = rowMtx->getVectorLengths()[r];
            for (int ki3 = 0; ki3 < rl3; ki3++) {
              int cc = rowMtx->getIndices()[rs3 + ki3];
              if (originalSolver->isInteger(cc))
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
          // Diagnostic: print details of remaining violations so we can
          // understand why the repair failed.
          fprintf(stderr, "[CbcRepairPostprocess] REPAIR FAILED: %d violation(s) remain.\n", bestViol);
          for (int r = 0; r < nRows; r++) {
            if (!pureInt[r])
              continue;
            double viol = std::max(rowLb[r] - lhs[r], lhs[r] - rowUb[r]);
            if (viol <= 1e-6)
              continue;
            fprintf(stderr, "  Row %d: lhs=%.10g lb=%.10g ub=%.10g viol=%.10g\n",
              r, lhs[r], rowLb[r], rowUb[r], viol);
            const int rs2 = rowMtx->getVectorStarts()[r];
            const int rl2 = rowMtx->getVectorLengths()[r];
            for (int ki2 = 0; ki2 < rl2; ki2++) {
              int c = rowMtx->getIndices()[rs2 + ki2];
              double coeff = rowMtx->getElements()[rs2 + ki2];
              double ppVal = (c < (int)postProcessSol.size()) ? postProcessSol[c] : -999.0;
              double repVal = repSol[c];
              bool retained = (c < (int)isRetainedEarly.size()) ? isRetainedEarly[c] : false;
              fprintf(stderr,
                "    col %d coeff=%.6g repSol=%.6g postProcess=%.6g retained=%d"
                " lb=%.6g ub=%.6g\n",
                c, coeff, repVal, ppVal, (int)retained,
                origColLb[c], origColUb[c]);
            }
          }
        }
      }

      // Update changed[] for any variable that moved from Phase-2 state.
      for (int c = 0; c < nCols; c++) {
        if (fabs(repSol[c] - repSolP3[c]) > 0.5)
          changed[c] = true;
      }
    }
  }

  // Inject the final repSol back into saveSolver.
  // This is needed because repSol may have been corrected via the
  // column mapping (overriding wrong back-substituted values) even
  // when no greedy repair flips occurred.
  // Additionally, fix all integer variables at their repSol values so
  // that the subsequent LP solves do not overwrite them.
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
