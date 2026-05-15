// Copyright (C) 2024 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#ifndef CbcBoundPropagation_hpp
#define CbcBoundPropagation_hpp

#include "CbcConfig.h"

class OsiSolverInterface;
class CoinMessageHandler;

/*! \brief Bound propagation: singleton tightening + knapsack-based
 *         bound propagation.
 *
 * This class orchestrates a two-phase bound propagation pipeline that runs
 * before the initial LP relaxation solve:
 *
 *  1. **Singleton tightening** (via OsiSolverInterface::tightenBoundsFromSingletonRows)
 *     — detects rows with a single binary/integer variable and fixes it when
 *     the constraint forces a unique value. Very fast: O(rows + cols).
 *
 *  2. **CoinBoundPropagation** — knapsack-based bound propagation over all
 *     ≤-type binary rows. Can cascade: after each round, newly fixed variables
 *     tighten other rows, enabling further fixings.
 *
 * ### Usage
 * \code
 *   CbcBoundPropagation bp;
 *   bool feasible = bp.run(solver, handler, logLevel,
 *                          CbcBoundPropagation::MILPbt,
 *                          100, useElapsed, timeLimit, startTime);
 *   if (!feasible) {
 *     // instance proved infeasible
 *     int row = bp.infeasibleRow();   // -1 if unknown
 *     int col = bp.infeasibleCol();   // -1 if unknown
 *   }
 *   // tightenings already applied to solver in place
 *   std::cout << "Fixed " << bp.nFixed() << " variables\n";
 * \endcode
 */
class CBCLIB_EXPORT CbcBoundPropagation {
public:
  /*! \brief Aggression level for bound propagation.
   *
   *  - Off:       nothing is done (caller should not call run()).
   *  - Singletons: singleton rows only.
   *  - MILPbt:   singletons then up to maxRounds of CoinBoundPropagation.
   *  - Fixpoint:  singletons then CoinBoundPropagation until fixpoint
   *               (ignores maxRounds).
   */
  enum Level {
    Off = 0,
    Singletons,
    MILPbt,
    Fixpoint
  };

  /*! \brief Reason why the bound propagation loop stopped. */
  enum StopReason {
    NotRun = 0,
    ReachedFixpoint, ///< No new fixings in last round
    HitMaxRounds, ///< Reached the maxRounds limit
    HitTimeLimit, ///< Ran out of time budget
    InfeasibleDetected ///< Proved infeasible
  };

  CbcBoundPropagation();

  /*! \brief Run bound propagation.
   *
   * Tightenings are applied in-place to \p solver.
   *
   * \param solver     The (already-loaded) MIP solver interface.
   * \param handler    Message handler for logging (may be null for silence).
   * \param logLevel   CbcModel log level (0=silent, 1=normal, 2=verbose).
   * \param level      Aggression level (must not be Off).
   * \param maxRounds  Maximum MILPbt rounds (only used when level == MILPbt).
   * \param useElapsed True to measure wall-clock time; false for CPU time.
   * \param timeLimit  Budget in seconds; use a large value to disable.
   * \param startTime  Reference time (from CoinGetTimeOfDay or CoinCpuTime)
   *                   measured BEFORE calling run().
   *
   * \return true if the problem is (still) feasible, false if infeasibility
   *         was proved during bound propagation.
   */
  bool run(OsiSolverInterface *solver,
    CoinMessageHandler *handler,
    int logLevel,
    Level level,
    int maxRounds,
    bool useElapsed,
    double timeLimit,
    double startTime);

  /// Number of variables tightened (not fully fixed) by singleton step.
  int nSingletonTightened() const { return nSingletonTightened_; }

  /// Number of variables fully fixed by singleton step.
  int nSingletonFixed() const { return nSingletonFixed_; }

  /// Number of variables fixed by CoinBoundPropagation across all rounds.
  int nBoundPropFixed() const { return nBoundPropFixed_; }

  /// Total fixings from all phases.
  int nFixed() const { return nSingletonFixed_ + nBoundPropFixed_; }

  /// Number of CoinBoundPropagation rounds executed.
  int nRoundsRun() const { return nRoundsRun_; }

  /// Reason the loop stopped (NotRun if run() was never called).
  StopReason stopReason() const { return stopReason_; }

  /// Elapsed/CPU time used by run() in seconds.
  double timeUsed() const { return timeUsed_; }

  /// Row index that caused infeasibility, or -1 if not available.
  int infeasibleRow() const { return infeasibleRow_; }

  /// Column index that caused infeasibility, or -1 if not applicable.
  int infeasibleCol() const { return infeasibleCol_; }

private:
  int nSingletonTightened_;
  int nSingletonFixed_;
  int nBoundPropFixed_;
  int nRoundsRun_;
  StopReason stopReason_;
  double timeUsed_;
  int infeasibleRow_;
  int infeasibleCol_;
};

#endif // CbcBoundPropagation_hpp
