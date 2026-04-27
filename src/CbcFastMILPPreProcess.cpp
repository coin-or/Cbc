// Copyright (C) 2024 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcFastMILPPreProcess.hpp"

#include "CoinMILPBoundTightening.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinTime.hpp"
#include "OsiSolverInterface.hpp"

#include <cassert>
#include <climits>
#include <cmath>
#include <string>
#include <vector>

CbcFastMILPPreProcess::CbcFastMILPPreProcess()
  : nSingletonTightened_(0)
  , nSingletonFixed_(0)
  , nMILPbtFixed_(0)
  , nRoundsRun_(0)
  , stopReason_(NotRun)
  , timeUsed_(0.0)
  , infeasibleRow_(-1)
  , infeasibleCol_(-1)
{
}

bool CbcFastMILPPreProcess::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel,
  Level level,
  int maxRounds,
  bool useElapsed,
  double timeLimit,
  double startTime)
{
  assert(level != Off);

  const double t0 = useElapsed ? CoinGetTimeOfDay() : CoinCpuTime();

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

      if (logLevel >= 1)
        printf("  Fast preprocessing: INFEASIBLE (singleton tightening), "
               "%.3f s.\n",
          timeUsed_);

      return false;
    }

    nSingletonFixed_ = nFixed;
    nSingletonTightened_ = nTightened - nFixed;

    if (logLevel >= 2 && nTightened > 0)
      printf("  Fast preprocessing (singletons): %d tightened, %d fixed.\n",
        nSingletonTightened_, nSingletonFixed_);
  }

  if (level == Singletons) {
    stopReason_ = ReachedFixpoint;
    timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

    if (logLevel >= 1)
      printf("  Fast preprocessing: fixed %d, tightened %d (singletons only), "
             "%.3f s.\n",
        nSingletonFixed_, nSingletonTightened_, timeUsed_);

    return true;
  }

  // ---------------------------------------------------------------
  // Phase 2: CoinMILPBoundTightening — iterate until fixpoint or limits
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
    const double tNow = useElapsed ? CoinGetTimeOfDay() : CoinCpuTime();
    if (tNow - startTime >= timeLimit) {
      stopReason_ = HitTimeLimit;
      timeUsed_ = tNow - t0;

      if (logLevel >= 1)
        printf("  Fast preprocessing: fixed %d (%d singleton + %d MILPbt, "
               "%d round(s), TIME LIMIT), %.3f s.\n",
          nSingletonFixed_ + nMILPbtFixed_, nSingletonFixed_, nMILPbtFixed_,
          nRoundsRun_, timeUsed_);

      return true;
    }

    // Construction runs the algorithm; results are immediately available.
    CoinMILPBoundTightening bt(nCols, colType,
      curLB.data(), curUB.data(),
      matByRow, rowSense, rhs, range,
      primalTol, infinity);

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
          printf("  Fast preprocessing: INFEASIBLE in round %d — "
                 "row %d (%s), col %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(),
            infeasibleCol_, colName.c_str(), timeUsed_);
        } else if (infeasibleRow_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          printf("  Fast preprocessing: INFEASIBLE in round %d — "
                 "row %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(), timeUsed_);
        } else {
          printf("  Fast preprocessing: INFEASIBLE in round %d, %.3f s.\n",
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
      timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

      if (logLevel >= 2)
        printf("  Fast preprocessing (MILPbt): fixpoint reached after %d "
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
    }

    nMILPbtFixed_ += nNew;

    if (logLevel >= 2)
      printf("  Fast preprocessing (MILPbt): round %d fixed %d variables "
             "(total %d).\n",
        nRoundsRun_, nNew, nMILPbtFixed_);
  }

  if (stopReason_ == NotRun) {
    // Exited by roundLimit without fixpoint
    stopReason_ = HitMaxRounds;
  }
  timeUsed_ = (useElapsed ? CoinGetTimeOfDay() : CoinCpuTime()) - t0;

  if (logLevel >= 1) {
    const int totalFixed = nSingletonFixed_ + nMILPbtFixed_;
    const char *stopNote = (stopReason_ == HitMaxRounds) ? ", MAX ROUNDS" : "";
    printf("  Fast preprocessing: fixed %d (%d singleton + %d MILPbt, "
           "%d round(s)%s), %.3f s.\n",
      totalFixed, nSingletonFixed_, nMILPbtFixed_,
      nRoundsRun_, stopNote, timeUsed_);
  }

  return true;
}
