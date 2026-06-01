// Copyright (C) 2024 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcBoundPropagation.hpp"

#include "CoinBoundPropagation.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"

#include <cassert>
#include <climits>
#include <cmath>
#include <string>
#include <vector>

CbcBoundPropagation::CbcBoundPropagation()
  : nSingletonTightened_(0)
  , nSingletonFixed_(0)
  , nBoundPropFixed_(0)
  , nRoundsRun_(0)
  , stopReason_(NotRun)
  , timeUsed_(0.0)
  , infeasibleRow_(-1)
  , infeasibleCol_(-1)
{
}

bool CbcBoundPropagation::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel,
  Level level,
  int maxRounds,
  double timeLimit,
  double startTime)
{
  assert(level != Off);

  if (!solver || solver->getNumCols() == 0)
    return true;

  const double t0 = CoinGetTimeOfDay();

  // If a debugCuts reference solution is loaded AND we are currently on the
  // optimal path (getRowCutDebugger returns non-NULL), check every bound
  // fixing we apply against that solution.  getRowCutDebugger (without
  // "Always") returns NULL at subtree nodes whose branching decisions have
  // already excluded the reference solution — correct behaviour there is to
  // fix variables differently from the global optimal, so we must NOT flag
  // those.  Only on the optimal path is a contradictory fixing a bug.
  const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
  const double *optSol = debugger ? debugger->optimalSolution() : nullptr;

  auto checkFixing = [&](int col, double lb, double ub,
                         const char *phase) {
    if (!optSol || !solver->isInteger(col))
      return;
    const double sv = optSol[col];
    if (lb > sv + 0.5 || ub < sv - 0.5) {
      const std::string name = solver->getColName(col);
      printf("nodeBoundProp BAD FIXING (%s): col %d (%s)"
             " fixed to [%g, %g] but optimal solution has %g\n",
             phase, col, name.c_str(), lb, ub, sv);
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
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 1)
        printf("  Bound propagation: INFEASIBLE (singleton tightening), "
               "%.3f s.\n",
          timeUsed_);

      return false;
    }

    nSingletonFixed_ = nFixed;
    nSingletonTightened_ = nTightened - nFixed;

    // Check singleton fixings against the reference solution (optimal path only).
    if (optSol) {
      const double *lb = solver->getColLower();
      const double *ub = solver->getColUpper();
      for (int j = 0; j < solver->getNumCols(); j++)
        checkFixing(j, lb[j], ub[j], "singleton");
    }

    if (logLevel >= 2 && nTightened > 0)
      printf("  Bound propagation (singletons): %d tightened, %d fixed.\n",
        nSingletonTightened_, nSingletonFixed_);
  }

  if (level == Singletons) {
    stopReason_ = ReachedFixpoint;
    timeUsed_ = (CoinGetTimeOfDay()) - t0;

    if (logLevel >= 1)
      printf("  Bound propagation fixed %d vars in %.3f s.\n",
        nSingletonFixed_, timeUsed_);

    return true;
  }

  // ---------------------------------------------------------------
  // Phase 2: CoinBoundPropagation — iterate until fixpoint or limits
  // ---------------------------------------------------------------
  const int roundLimit = (level == Fixpoint) ? INT_MAX : maxRounds;

  // Extract solver data once; bounds are updated in the loop.
  const int nCols = solver->getNumCols();
  const char *colType = solver->getColType(false);
  const CoinPackedMatrix *matByRow = solver->getMatrixByRow();
  const char *rowSense = solver->getRowSense();
  const double *rhs = solver->getRightHandSide();
  const double *range = solver->getRowRange();
  double primalTol = 1e-7;
  solver->getDblParam(OsiPrimalTolerance, primalTol);
  const double infinity = solver->getInfinity();

  // Working copies of bounds (updated after each round).
  std::vector< double > curLB(solver->getColLower(),
    solver->getColLower() + nCols);
  std::vector< double > curUB(solver->getColUpper(),
    solver->getColUpper() + nCols);

  for (int round = 0; round < roundLimit; ++round) {
    // Time-limit check before starting this round
    const double tNow = CoinGetTimeOfDay();
    if (tNow - startTime >= timeLimit) {
      stopReason_ = HitTimeLimit;
      timeUsed_ = tNow - t0;

      if (logLevel >= 1)
        printf("  Bound propagation: fixed %d (%d singleton + %d propagation, "
               "%d round(s), TIME LIMIT), %.3f s.\n",
          nSingletonFixed_ + nBoundPropFixed_, nSingletonFixed_, nBoundPropFixed_,
          nRoundsRun_, timeUsed_);

      return true;
    }

    // Construction runs the algorithm; results are immediately available.
    CoinBoundPropagation bt(nCols, colType,
      curLB.data(), curUB.data(),
      matByRow, rowSense, rhs, range,
      primalTol, infinity);

    ++nRoundsRun_;

    if (bt.isInfeasible()) {
      infeasibleRow_ = bt.infeasibleRow();
      infeasibleCol_ = bt.infeasibleCol();
      stopReason_ = InfeasibleDetected;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 1) {
        if (infeasibleRow_ >= 0 && infeasibleCol_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          const std::string colName = (infeasibleCol_ < solver->getNumCols())
            ? solver->getColName(infeasibleCol_)
            : "(unknown)";
          printf("  Bound propagation: INFEASIBLE in round %d — "
                 "row %d (%s), col %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(),
            infeasibleCol_, colName.c_str(), timeUsed_);
        } else if (infeasibleRow_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          printf("  Bound propagation: INFEASIBLE in round %d — "
                 "row %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(), timeUsed_);
        } else {
          printf("  Bound propagation: INFEASIBLE in round %d, %.3f s.\n",
            nRoundsRun_, timeUsed_);
        }
      }

      return false;
    }

    // Apply the fixings from this round to curLB/curUB and to the solver
    const auto &bounds = bt.updatedBounds();
    const int nNew = static_cast< int >(bounds.size());

    if (nNew == 0) {
      stopReason_ = ReachedFixpoint;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 2)
        printf("  Bound propagation: fixpoint reached after %d "
               "round(s).\n",
          nRoundsRun_);

      break;
    }

    for (const auto &p : bounds) {
      const int col = static_cast< int >(p.first);
      curLB[col] = p.second.first;
      curUB[col] = p.second.second;
      solver->setColLower(col, p.second.first);
      solver->setColUpper(col, p.second.second);
      checkFixing(col, p.second.first, p.second.second, "propagation");
    }

    nBoundPropFixed_ += nNew;

    if (logLevel >= 2)
      printf("  Bound propagation: round %d fixed %d variables "
             "(total %d).\n",
        nRoundsRun_, nNew, nBoundPropFixed_);
  }

  if (stopReason_ == NotRun) {
    // Exited by roundLimit without fixpoint
    stopReason_ = HitMaxRounds;
  }
  timeUsed_ = (CoinGetTimeOfDay()) - t0;

  if (logLevel >= 1) {
    const int totalFixed = nSingletonFixed_ + nBoundPropFixed_;
    printf("  Bound propagation fixed %d vars in %.3f s.\n",
      totalFixed, timeUsed_);
  }

  return true;
}
