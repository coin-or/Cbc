// Copyright (C) 2024 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcBoundPropagation.hpp"

#include "CoinBoundPropagation.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"

// Reference solution loaded by -debugCuts; set before applyLpMethod() so that
// bound propagation can check fixings before the OsiRowCutDebugger is active.
extern double *debugSolution;
extern int debugNumberColumns;

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <string>
#include <vector>

CbcBoundPropagation::CbcBoundPropagation()
  : nSingletonTightened_(0)
  , nSingletonFixed_(0)
  , nBoundPropFixed_(0)
  , nFBBTTightened_(0)
  , nRoundsRun_(0)
  , stopReason_(NotRun)
  , timeUsed_(0.0)
  , infeasibleRow_(-1)
  , infeasibleCol_(-1)
  , nonBinaryFBBT_(true)
{
}

bool CbcBoundPropagation::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel,
  Level level,
  int maxRounds,
  bool useElapsed,
  double timeLimit,
  double startTime)
{
  assert(level != Off);

  if (!solver || solver->getNumCols() == 0)
    return true;

  const double t0 = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime());

  // If a debugCuts reference solution is loaded, check every bound fixing we
  // apply against that solution.
  //
  // Before the LP solve the OsiRowCutDebugger is not yet active, so we fall
  // back to the debugSolution global which is populated from the -debugCuts
  // sol file at the very start of applyLpMethod().
  //
  // NOTE: previously this used getRowCutDebugger() (without "Always"), which
  // relies on OsiRowCutDebugger::onOptimalPath() -- that check only looks at
  // INTEGER columns with a loose 1e-3 tolerance, and has been observed to be
  // unreliable (can return false even when the current bounds are still
  // fully consistent with the reference solution, and vice versa). Use
  // getRowCutDebuggerAlways() instead (bypasses onOptimalPath() entirely)
  // and compute our own manual on-path check across ALL columns with a
  // tight 1e-6 tolerance, matching the approach used elsewhere in
  // CbcModel.cpp (pre-resolve cut/column-cut checks, reducedCostFix(),
  // fixFromGlobalCuts()) for exactly this reason.
  const int nCols = solver->getNumCols();
  const OsiRowCutDebugger *debuggerAlways = solver->getRowCutDebuggerAlways();
  const double *optSol = nullptr;
  if (debuggerAlways) {
    const double *refSol = debuggerAlways->optimalSolution();
    const double *nowLB = solver->getColLower();
    const double *nowUB = solver->getColUpper();
    bool stillOnPath = true;
    for (int k = 0; k < nCols; k++) {
      if (refSol[k] < nowLB[k] - 1.0e-6 || refSol[k] > nowUB[k] + 1.0e-6) {
        stillOnPath = false;
        break;
      }
    }
    if (stillOnPath)
      optSol = refSol;
  } else if (debugSolution && debugNumberColumns == solver->getNumCols()) {
    // Pre-B&B bound propagation: no debugger activated yet, fall back to the
    // raw -debugCuts reference solution global.
    optSol = debugSolution;
  }

  // Declare these early so they are in scope for the checkFixing lambda.
  // colType and curLB/curUB are set to their real values before phase 2.
  // Local owned copy — not a pointer into the solver's internal cache.
  // Refreshed before each propagation round from current bounds.
  std::vector< char > colTypeBuf;
  const char *colType = nullptr;       // points to colTypeBuf.data() after phase-2 init
  std::vector< double > curLB, curUB;  // initialised before phase-2 loop

  // diagRound == -1: singleton phase;  >= 0: CoinBoundPropagation round.
  int diagRound = -1;

  // Check one fixing against the reference solution (optimal path only).
  // Called with the new bounds; curLB/curUB must still hold OLD bounds when
  // called from the propagation loop (call before updating curLB/curUB).
  auto checkFixing = [&](int col, double newLB, double newUB,
                         const char *phase) {
    if (!optSol)
      return;
    const double sv = optSol[col];
    const bool isInt = solver->isInteger(col);
    // For integers: use 0.5 tolerance (any integer rounding is wrong).
    // For continuous: use 1e-8 tolerance (catches FBBT bounds tighter than
    //   the true mathematical bound, which exclude the reference solution).
    const double tol = isInt ? 0.5 : 1.0e-8;
    if (newLB > sv + tol || newUB < sv - tol) {
      const std::string name = solver->getColName(col);
      const char ct = colType ? colType[col] : char(-1);
      const char *ctName = (ct == 1) ? "binary"
                         : (ct == 2) ? "general-integer"
                         : (ct == 0) ? "continuous" : "?";
      // curLB/curUB contain old bounds when called from propagation loop.
      const double oldLB = (col < static_cast< int >(curLB.size()))
                             ? curLB[col] : newLB;
      const double oldUB = (col < static_cast< int >(curUB.size()))
                             ? curUB[col] : newUB;
      printf("nodeBoundProp BAD FIXING (%s, round %d): col %d (%s)"
             " type=%s old=[%g,%g] new=[%g,%g] but optimal has %g\n",
             phase, diagRound, col, name.c_str(), ctName,
             oldLB, oldUB, newLB, newUB, sv);
      fflush(stdout);
    }
  };

  // ---------------------------------------------------------------
  // Phase 1: singleton row tightening
  // ---------------------------------------------------------------
  {
    int nFixed = 0;
    const int nTightened = solver->tightenBoundsFromSingletonRows(nFixed);

    if (nTightened < 0) {
      // infeasibility detected — singleton API does not expose row/col source
      stopReason_ = InfeasibleDetected;
      timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

      if (logLevel >= 1) {
        printf("  Bound tightening: INFEASIBLE (singleton rows), "
               "%.3f s.\n",
          timeUsed_);
        fflush(stdout);
      }

      return false;
    }

    nSingletonFixed_ = nFixed;
    nSingletonTightened_ = nTightened - nFixed;

    // Check singleton fixings against the reference solution (optimal path only).
    // curLB/curUB are not yet populated; oldLB/oldUB will show newLB/newUB.
    if (optSol) {
      // Capture stable copies — the raw pointers from getColLower/Upper() may
      // be invalidated by intermediate solver calls (e.g. getColName(),
      // isInteger()) that trigger internal array reallocation.
      const double *rawLB = solver->getColLower();
      const double *rawUB = solver->getColUpper();
      const std::vector<double> lbVec(rawLB, rawLB + nCols);
      const std::vector<double> ubVec(rawUB, rawUB + nCols);
      for (int j = 0; j < nCols; j++)
        checkFixing(j, lbVec[j], ubVec[j], "singleton");
    }

    if (logLevel >= 2 && nTightened > 0) {
      printf("  Bound propagation (singletons): %d tightened, %d fixed.\n",
        nSingletonTightened_, nSingletonFixed_);
      fflush(stdout);
    }
  }

  if (level == Singletons) {
    stopReason_ = ReachedFixpoint;
    timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

    // Singleton tightening and fixpoint/FBBT propagation are all just
    // row-based bound tightening — report a single unified count instead of
    // separate "bound propagation"/"FBBT" terminology (a fixed variable is
    // simply the extreme case of a tightened bound, so it is included in the
    // total).
    if (logLevel >= 1) {
      printf("  Bound tightening: %d bounds tightened (%d fixed) in %.3f s.\n",
        nSingletonTightened_ + nSingletonFixed_, nSingletonFixed_, timeUsed_);
      fflush(stdout);
    }

    return true;
  }

  // ---------------------------------------------------------------
  // Phase 2: CoinBoundPropagation — iterate until fixpoint or limits
  // ---------------------------------------------------------------

  // Initialise phase-2 data (after singleton tightening so bounds are current).
  // Copy colType into a local buffer so we own the data (not a pointer into
  // the solver's internal cache which may be stale from a prior invocation).
  // Use refresh=true to force recomputation from current bounds — a stale
  // cached value can misclassify GeneralInteger variables as Binary, causing
  // wrong fixings (root cause of the miclsp arm64 soundness bug).
  {
    const char *ct = solver->getColType(true);
    colTypeBuf.assign(ct, ct + nCols);
  }
  colType = colTypeBuf.data();
  const CoinPackedMatrix *matByRow = solver->getMatrixByRow();
  const char *rowSense = solver->getRowSense();
  const double *rhs = solver->getRightHandSide();
  const double *range = solver->getRowRange();
  double primalTol = 1e-7;
  solver->getDblParam(OsiPrimalTolerance, primalTol);
  const double infinity = solver->getInfinity();

  // Working copies of bounds (updated after each round).
  curLB.assign(solver->getColLower(), solver->getColLower() + nCols);
  curUB.assign(solver->getColUpper(), solver->getColUpper() + nCols);

  // Refresh colTypeBuf from current bounds before each round.
  // Cheaper than calling solver->getColType() and avoids stale-cache issues.
  // A GeneralInteger whose bounds have been tightened to [0,1] is reclassified
  // as Binary, enabling stronger propagation in later rounds.
  auto refreshColType = [&]() {
    for (int j = 0; j < nCols; ++j) {
      if (colTypeBuf[j] == 0)
        continue; // continuous: never changes
      colTypeBuf[j] = (curLB[j] >= 0.0 && curUB[j] <= 1.0) ? 1 : 2;
    }
  };

  // ── Dirty-row tracking for FBBT (lazy) ───────────────────────────────────
  // Round 1 always runs without a dirty hint (CoinBoundPropagation builds its
  // own rowHasNonBinary).  The O(nnz) colToRows adjacency list and the dirty
  // bitvector are built lazily after round 1 — only if there will be a round 2.
  // This avoids paying the build cost for problems that converge in one round.
  const int nRows = matByRow->getNumRows();
  const int *matIdxs = matByRow->getIndices();
  const CoinBigIndex *matStart = matByRow->getVectorStarts();
  const int *matLen = matByRow->getVectorLengths();

  bool hasNonBinaryVars = false;
  if (nonBinaryFBBT_) {
    for (int j = 0; j < nCols && !hasNonBinaryVars; ++j)
      if (colTypeBuf[j] != 1) // 1 = Binary
        hasNonBinaryVars = true;
  }

  // ---------------------------------------------------------------
  // Divergence (positive-cycle) detection for FBBT.
  // ---------------------------------------------------------------
  // Rows of the form  x_j - x_i >= c  (disjunctive/precedence constraints,
  // e.g. job-shop scheduling) implement, under FBBT, the exact same bound
  // relaxation rule as Bellman-Ford shortest-path relaxation:
  //   lb(x_j) = max(lb(x_j), lb(x_i) + c).
  // Bellman-Ford's classic guarantee is that shortest paths (here: tightest
  // finite lower bounds) stabilise within at most (nCols - 1) full relaxation
  // passes UNLESS the underlying graph contains a positive-weight cycle, in
  // which case relaxation improves some bound forever. That is a genuine
  // MATHEMATICAL CERTIFICATE OF INFEASIBILITY for a difference-constraint
  // subsystem — a positive cycle means the constraints are contradictory
  // (e.g. x1 - x0 >= 5 and x0 - x1 >= 5 imply 0 >= 10) — NOT merely "FBBT
  // gave up too early". Our round-robin, dirty-row FBBT sweep is a Jacobi-style
  // generalisation of the same relaxation (rows may touch more than 2
  // variables), and the same argument applies: if round-to-round tightening
  // has not decayed after touring the full variable set once (nCols rounds)
  // it cannot be legitimate slow convergence — a converging system must show
  // shrinking deltas well within that many rounds — so it must be a cycle.
  //
  // We detect this directly (not by an arbitrary hard round cap): track the
  // total tightening magnitude (sum of |ΔLB|+|ΔUB|) applied per round. Once
  // more than nCols rounds have elapsed with no fixings (a fixing legitimately
  // changes the system and resets the window) and the magnitude has not
  // decayed by at least half compared to nCols rounds earlier, declare
  // infeasibility right away — this is both CORRECT (matches what the LP
  // resolve that follows would eventually discover via a Farkas ray) and much
  // FASTER (fires within ~nCols rounds instead of spinning indefinitely).
  std::vector< double > roundDeltaHistory;
  const int divergenceWindow = std::max(10, nCols);

  // Fixpoint level has no explicit maxRounds cap by design (see header doc:
  // "ignores maxRounds") since it targets the genuine mathematical fixpoint.
  // The divergence detector above is expected to terminate any non-converging
  // case within O(nCols) rounds; roundLimit itself is only a last-resort
  // safety net (never expected to be hit) in case some pathological instance
  // manages to evade the detector.
  const int roundLimit = (level == Fixpoint) ? INT_MAX : maxRounds;

  // Built lazily after round 1 (if a second round is needed).
  bool dirtyInfraBuilt = false;
  std::vector< CoinBigIndex > colRowStart;
  std::vector< int > colRowList;
  std::vector< double > colRowCoef; // coefficient of col j in each (col,row) entry
  std::vector< bool > rowHasNonBinaryBP;
  std::vector< char > rowHasBinaryBP; // char to allow .data() (vector<bool> lacks it)
  std::vector< bool > dirtyRowsFBBT;

  // Cached row min activity (for ≤ constraint FBBT, incremental updates).
  // rowCachedMinAct[r] = sum(a>0: a*lb_j, a<0: a*ub_j) over all vars in row r.
  // rowCachedNUnbLB[r] = number of vars with a*lb_j = -inf (contributes -inf to minAct).
  // Valid when actCacheBuilt == true; updated incrementally from bound changes.
  bool actCacheBuilt = false;
  std::vector< double > rowCachedMinAct;
  std::vector< int > rowCachedNUnbLB;

  // Per-round buffer recording (col, oldLB, oldUB, newLB, newUB) for each committed
  // bound change. Used in the post-round single adjacency pass that simultaneously
  // updates the min-activity cache and marks dirty rows. Pre-allocated once, reused
  // every round to avoid repeated heap allocations.
  struct BoundChange {
    int col;
    double oldLB, oldUB, newLB, newUB;
  };
  std::vector< BoundChange > changedBounds;

  // Build colToRows (column → rows adjacency) and rowHasNonBinaryBP once.
  // Called at most once, the first time a round 2 is needed.
  auto buildDirtyInfra = [&]() {
    if (dirtyInfraBuilt)
      return;
    dirtyInfraBuilt = true;

    rowHasNonBinaryBP.assign(static_cast< size_t >(nRows), false);
    rowHasBinaryBP.assign(static_cast< size_t >(nRows), char(0));
    for (int r = 0; r < nRows; ++r) {
      const CoinBigIndex rs = matStart[r];
      const int len = matLen[r];
      for (int k = 0; k < len; ++k) {
        const int j = matIdxs[rs + k];
        if (colTypeBuf[j] != 1)
          rowHasNonBinaryBP[r] = true;
        else
          rowHasBinaryBP[r] = char(1);
      }
    }

    // colToRows: for each column j, list the rows and coefficients.
    colRowStart.assign(static_cast< size_t >(nCols + 1), CoinBigIndex(0));
    for (int r = 0; r < nRows; ++r) {
      const CoinBigIndex rs = matStart[r];
      const int len = matLen[r];
      for (int k = 0; k < len; ++k)
        ++colRowStart[matIdxs[rs + k] + 1];
    }
    for (int j = 0; j < nCols; ++j)
      colRowStart[j + 1] += colRowStart[j];
    colRowList.resize(static_cast< size_t >(colRowStart[nCols]));
    colRowCoef.resize(static_cast< size_t >(colRowStart[nCols]));
    {
      std::vector< CoinBigIndex > pos(colRowStart.begin(),
        colRowStart.begin() + nCols);
      const double *matCoefs = matByRow->getElements();
      for (int r = 0; r < nRows; ++r) {
        const CoinBigIndex rs = matStart[r];
        const int len = matLen[r];
        for (int k = 0; k < len; ++k) {
          const int j = matIdxs[rs + k];
          const CoinBigIndex p = pos[j]++;
          colRowList[p] = r;
          colRowCoef[p] = matCoefs[rs + k];
        }
      }
    }

    // Initialise the row min-activity cache from the current bounds (curLB/curUB).
    // Use Kahan compensated summation to match the precision of scanRow().
    // Without this, plain sequential FP accumulation on arm64 at -O1 can compute
    // minAct slightly above the true value, causing the same over-tight FBBT
    // bounds as in scanRow (the cache bypasses scanRow from round 2 onwards).
    rowCachedMinAct.assign(static_cast< size_t >(nRows), 0.0);
    rowCachedNUnbLB.assign(static_cast< size_t >(nRows), 0);
    {
      const double *matCoefs = matByRow->getElements();
      for (int r = 0; r < nRows; ++r) {
        const CoinBigIndex rs = matStart[r];
        const int len = matLen[r];
        double minA = 0.0, minAComp = 0.0;
        int nUnbLB = 0;
        for (int k = 0; k < len; ++k) {
          const int j = matIdxs[rs + k];
          const double a = matCoefs[rs + k];
          const double lb = curLB[j], ub = curUB[j];
          double term = 0.0;
          if (a > 0.0) {
            if (lb <= -infinity)
              ++nUnbLB;
            else
              term = a * lb;
          } else if (a < 0.0) {
            if (ub >= infinity)
              ++nUnbLB;
            else
              term = a * ub;
          }
          // Kahan compensated addition
          const double y = term - minAComp;
          const double t = minA + y;
          minAComp = (t - minA) - y;
          minA = t;
        }
        rowCachedMinAct[r] = minA;
        rowCachedNUnbLB[r] = nUnbLB;
      }
    }
    actCacheBuilt = true;

    dirtyRowsFBBT.assign(static_cast< size_t >(nRows), false);
  };

  for (int round = 0; round < roundLimit; ++round) {
    diagRound = round;
    refreshColType(); // update binary/general-integer from current curLB/curUB
    // Time-limit check before starting this round
    const double tNow = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime());
    if (tNow - startTime >= timeLimit) {
      stopReason_ = HitTimeLimit;
      timeUsed_ = tNow - t0;

      if (logLevel >= 1) {
        printf("  Bound tightening: fixed %d (%d singleton + %d propagation, "
               "%d round(s), TIME LIMIT), %.3f s.\n",
          nSingletonFixed_ + nBoundPropFixed_, nSingletonFixed_, nBoundPropFixed_,
          nRoundsRun_, timeUsed_);
        fflush(stdout);
      }

      return true;
    }

    // Round 0 and 1 always use nullptr (no dirty hint; CoinBoundPropagation
    // builds rowHasNonBinary internally).  Round 2+ use the dirty set built
    // after round 1 completes.
    const std::vector< bool > *dirtyHint =
      (hasNonBinaryVars && dirtyInfraBuilt) ? &dirtyRowsFBBT : nullptr;
    const double *cachedMinAct = actCacheBuilt ? rowCachedMinAct.data() : nullptr;
    const int *cachedNUnbLB = actCacheBuilt ? rowCachedNUnbLB.data() : nullptr;
    const bool *hasBinaryRow =
      dirtyInfraBuilt ? reinterpret_cast< const bool * >(rowHasBinaryBP.data()) : nullptr;
    // Heartbeat before the (potentially expensive, non-preemptible) scan
    // below — this call has no internal time-limit check, so on a large or
    // pathological instance it can run well past the deadline in a single
    // shot. Printing+flushing here lets a killed run's log pinpoint whether
    // it died mid-round (this line is the last one) or elsewhere.
    if (logLevel >= 2) {
      printf("  Bound propagation: round %d starting (nRows=%d nCols=%d "
             "dirtyHint=%s), t=%.3fs\n",
        round, nRows, nCols, dirtyHint ? "yes" : "no", tNow - startTime);
      fflush(stdout);
    }
    CoinBoundPropagation bt(nCols, colType,
      curLB.data(), curUB.data(),
      matByRow, rowSense, rhs, range,
      primalTol, infinity,
      /*maxRowNz=*/-1, /*collectCases=*/false,
      nonBinaryFBBT_, dirtyHint,
      cachedMinAct, nullptr, cachedNUnbLB, nullptr,
      hasBinaryRow);

    ++nRoundsRun_;

    if (bt.isInfeasible()) {
      infeasibleRow_ = bt.infeasibleRow();
      infeasibleCol_ = bt.infeasibleCol();
      stopReason_ = InfeasibleDetected;
      timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

      if (logLevel >= 1) {
        if (infeasibleRow_ >= 0 && infeasibleCol_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          const std::string colName = (infeasibleCol_ < solver->getNumCols())
            ? solver->getColName(infeasibleCol_)
            : "(unknown)";
          printf("  Bound tightening: INFEASIBLE in round %d — "
                 "row %d (%s), col %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(),
            infeasibleCol_, colName.c_str(), timeUsed_);
        } else if (infeasibleRow_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          printf("  Bound tightening: INFEASIBLE in round %d — "
                 "row %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(), timeUsed_);
        } else {
          printf("  Bound tightening: INFEASIBLE in round %d, %.3f s.\n",
            nRoundsRun_, timeUsed_);
        }
        fflush(stdout);
      }

      return false;
    }

    // Apply the fixings from this round to curLB/curUB and to the solver.
    // newBounds_ contains binary fixings (first) then FBBT tightenings (appended
    // at end of constructor), but we apply all in one loop.
    const auto &bounds = bt.updatedBounds();
    const int nNew = static_cast< int >(bounds.size());
    const int nFBBT = bt.nContinuousTightened();
    const int nFixed = nNew - nFBBT;

    if (nNew == 0) {
      stopReason_ = ReachedFixpoint;
      timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

      if (logLevel >= 2) {
        printf("  Bound propagation: fixpoint reached after %d "
               "round(s).\n",
          nRoundsRun_);
        fflush(stdout);
      }

      break;
    }

    // Commit bound changes and record old/new for the post-round adjacency pass.
    changedBounds.clear();
    double sumAbsDeltaThisRound = 0.0;
    int worstCol = -1;
    double worstColDelta = 0.0;
    for (const auto &p : bounds) {
      const int col = static_cast< int >(p.first);
      const double newLB = p.second.first;
      const double newUB = p.second.second;
      const double oldLB = curLB[col];
      const double oldUB = curUB[col];
      // Check BEFORE updating curLB/curUB so the lambda sees the old bounds.
      checkFixing(col, newLB, newUB, "propagation");
      // logLevel >= 3: trace every individual bound change with its delta, so
      // a genuine (but very slow, geometrically-converging) tightening chain
      // can be told apart from a bug that keeps re-touching the same column
      // by a non-shrinking or oscillating amount round after round.
      if (logLevel >= 3 && (newLB != oldLB || newUB != oldUB)) {
        printf("  [BP diag] round %d: col %d dLB=%.3e (%.15g->%.15g)"
               " dUB=%.3e (%.15g->%.15g)\n",
          round, col, newLB - oldLB, oldLB, newLB, newUB - oldUB, oldUB, newUB);
        fflush(stdout);
      }
      // Only accumulate finite-to-finite deltas — a bound moving between two
      // "infinity" sentinels is not a real tightening.
      const double absDLB = (oldLB > -infinity && newLB > -infinity) ? std::fabs(newLB - oldLB) : 0.0;
      const double absDUB = (oldUB < infinity && newUB < infinity) ? std::fabs(newUB - oldUB) : 0.0;
      const double colDelta = absDLB + absDUB;
      sumAbsDeltaThisRound += colDelta;
      if (colDelta > worstColDelta) {
        worstColDelta = colDelta;
        worstCol = col;
      }
      curLB[col] = newLB;
      curUB[col] = newUB;
      solver->setColLower(col, newLB);
      solver->setColUpper(col, newUB);
      if (hasNonBinaryVars)
        changedBounds.push_back({ col, oldLB, oldUB, newLB, newUB });
    }

    nBoundPropFixed_ += nFixed;
    nFBBTTightened_ += nFBBT;

    // Post-round single adjacency pass: update min-activity cache AND mark dirty rows.
    // Deferred until round >= 1 so fast 1-2 round problems never pay the build cost.
    if (hasNonBinaryVars && round >= 1) {
      const bool cacheAlreadyBuilt = actCacheBuilt; // false on first call (round 1)
      buildDirtyInfra(); // builds cache from current curLB/curUB (first time only)
      std::fill(dirtyRowsFBBT.begin(), dirtyRowsFBBT.end(), false);
      for (const BoundChange &bc : changedBounds) {
        const double dLB = bc.newLB - bc.oldLB;
        const double dUB = bc.newUB - bc.oldUB;
        for (CoinBigIndex ri = colRowStart[bc.col];
             ri < colRowStart[bc.col + 1]; ++ri) {
          const int r = colRowList[ri];
          // Update min-activity cache — only when the cache was already built
          // BEFORE this round's changes (i.e., not on the initial build round).
          // On the initial build, curLB/curUB already include round 1's changes,
          // so applying the deltas again would double-count.
          if (cacheAlreadyBuilt) {
            const double a = colRowCoef[ri];
            if (a > 0.0 && dLB != 0.0) {
              if (bc.oldLB <= -infinity) {
                --rowCachedNUnbLB[r];
                rowCachedMinAct[r] += a * bc.newLB;
              } else {
                rowCachedMinAct[r] += a * dLB;
              }
            } else if (a < 0.0 && dUB != 0.0) {
              if (bc.oldUB >= infinity) {
                --rowCachedNUnbLB[r];
                rowCachedMinAct[r] += a * bc.newUB;
              } else {
                rowCachedMinAct[r] += a * dUB;
              }
            }
          }
          // Mark dirty rows for next round.
          if (rowHasNonBinaryBP[r])
            dirtyRowsFBBT[r] = true;
        }
      }
    } // if (hasNonBinaryVars && round >= 1)

    if (logLevel >= 2) {
      printf("  Bound propagation: round %d fixed %d vars, FBBT tightened %d"
             " (total fixed %d).\n",
        nRoundsRun_, nFixed, nFBBT, nBoundPropFixed_);
      fflush(stdout);
    }

    // Divergence (positive-cycle) check — see the detailed rationale in the
    // setup block above where roundDeltaHistory/divergenceWindow are defined.
    if (nFixed > 0) {
      // A genuine fixing changes the system meaningfully; any prior "no
      // decay" history is no longer a valid baseline for comparison.
      roundDeltaHistory.clear();
    } else {
      roundDeltaHistory.push_back(sumAbsDeltaThisRound);
      if (static_cast< int >(roundDeltaHistory.size()) > divergenceWindow) {
        const double deltaNow = roundDeltaHistory.back();
        const double deltaBefore = roundDeltaHistory[roundDeltaHistory.size() - 1 - divergenceWindow];
        // Require at least 50% decay over a full window; anything less is
        // treated as non-convergent (tolerates floating-point noise while
        // still catching true divergence, which empirically shows exactly
        // zero decay round after round).
        if (deltaBefore > primalTol && deltaNow >= 0.5 * deltaBefore) {
          // Best-effort witness row: any currently-dirty row incident to the
          // worst-drifting column.
          int witnessRow = -1;
          if (dirtyInfraBuilt && worstCol >= 0) {
            for (CoinBigIndex ri = colRowStart[worstCol];
                 ri < colRowStart[worstCol + 1]; ++ri) {
              const int r = colRowList[ri];
              if (dirtyRowsFBBT[r]) {
                witnessRow = r;
                break;
              }
            }
          }
          infeasibleCol_ = worstCol;
          infeasibleRow_ = witnessRow;
          stopReason_ = InfeasibleDetected;
          timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

          if (logLevel >= 1) {
            const std::string colName = (infeasibleCol_ >= 0 && infeasibleCol_ < solver->getNumCols())
              ? solver->getColName(infeasibleCol_)
              : "(unknown)";
            if (infeasibleRow_ >= 0) {
              const std::string rowName = (infeasibleRow_ < solver->getNumRows())
                ? solver->getRowName(infeasibleRow_)
                : "(unknown)";
              printf("  Bound tightening: INFEASIBLE -- FBBT divergence detected"
                     " after %d round(s) with no decaying progress"
                     " (col %d (%s) still moving by %.3g per round,"
                     " row %d (%s), %.3f s -- a Bellman-Ford"
                     " positive-cycle certificate).\n",
                static_cast< int >(roundDeltaHistory.size()), infeasibleCol_,
                colName.c_str(), deltaNow, infeasibleRow_, rowName.c_str(), timeUsed_);
            } else {
              printf("  Bound tightening: INFEASIBLE -- FBBT divergence detected"
                     " after %d round(s) with no decaying progress"
                     " (col %d (%s) still moving by %.3g per round), %.3f s"
                     " -- a Bellman-Ford positive-cycle certificate.\n",
                static_cast< int >(roundDeltaHistory.size()), infeasibleCol_,
                colName.c_str(), deltaNow, timeUsed_);
            }
            fflush(stdout);
          }
          return false;
        }
      }
    }
  }

  if (stopReason_ == NotRun) {
    // Exited by roundLimit without fixpoint
    stopReason_ = HitMaxRounds;
  }

  timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

  if (logLevel >= 1) {
    // Singleton tightening, fixpoint propagation, and FBBT are all just
    // row-based bound tightening — report a single unified count rather than
    // separate "bound propagation"/"FBBT" terminology (a fixed variable is
    // simply the extreme case of a tightened bound, so it is included here).
    const int totalFixed = nSingletonFixed_ + nBoundPropFixed_;
    const int totalTightened = nSingletonTightened_ + totalFixed + nFBBTTightened_;
    printf("  Bound tightening: %d bounds tightened (%d fixed) in %.3f s.\n",
      totalTightened, totalFixed, timeUsed_);
    fflush(stdout);
  }

  return true;
}
