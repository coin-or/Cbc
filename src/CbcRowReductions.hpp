// Copyright (C) 2026 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#ifndef CbcRowReductions_hpp
#define CbcRowReductions_hpp

#include "CbcConfig.h"

class OsiSolverInterface;
class CoinMessageHandler;

/*! \brief Fast row reductions for the pre-root-LP strengthening phase.
 *
 * Deletes rows that cannot constrain the problem, before the root LP is
 * solved. Three reductions, one scan and one sort between them:
 *
 *  1. **All-fixed rows.** Every column of the row is fixed, so the row's
 *     activity is a constant: it is either violated (the model is infeasible)
 *     or it can never be binding (the row goes).
 *  2. **Duplicate rows.** Row \f$i\f$ and row \f$j\f$ have identical support
 *     and identical coefficients.
 *  3. **Parallel (scaled) rows.** Row \f$i = \lambda\,\text{row}_j\f$ for some
 *     \f$\lambda \ne 0\f$. The surviving row's bounds become the intersection
 *     of its own and the other row's, scaled by \f$\lambda\f$.
 *
 * 2 and 3 are the *same* reduction — an exact duplicate is just
 * \f$\lambda = 1\f$ — and are implemented as one pass. That is also how HiGHS
 * (HPresolve::detectParallelRowsAndCols) and papilo (ParallelRowDetection) do
 * it; neither has a separate duplicate-row pass. The counters are reported
 * separately only because the distinction is interesting, not because the code
 * paths differ.
 *
 * ### Why this is cheap
 *
 * The naive form of 2/3 is \f$O(m^2)\f$ row comparisons. Instead each row is
 * reduced to a hash that is *invariant under scaling*: divide every coefficient
 * by the row's pivot (the largest-magnitude entry, ties broken by smallest
 * column index) before hashing it, so \f$\lambda\,\text{row}\f$ and
 * \f$\text{row}\f$ hash identically. Sorting the row indices by that hash
 * groups all candidates into contiguous buckets, and only rows inside a bucket
 * are ever compared: \f$O(\text{nnz})\f$ to hash, \f$O(m \log m)\f$ to sort.
 *
 * Two details of the hash are load-bearing rather than cosmetic:
 *
 *  - **The pivot tie-break must be deterministic** (smallest column index
 *    among entries of equal largest magnitude). Without it two rows related by
 *    \f$\lambda = -1\f$ can pick different pivots, normalise to different
 *    vectors, and never land in the same bucket. Note this is a guarantee about
 *    the *interface*, not about any one implementation of it: OsiClp happens to
 *    store every row in ascending column order, which makes the tie-break a
 *    no-op there, but CoinPackedMatrix promises no such order.
 *  - **Coefficients are bucketed, not hashed exactly.** A value is multiplied
 *    by \f$1/\varphi\f$ and only its exponent and leading mantissa bits are
 *    kept, so two coefficients that differ by rounding noise hash the same
 *    instead of straddling a power-of-two boundary. The irrational multiplier
 *    is what keeps 0.5 and 0.5-\f$\epsilon\f$ on the same side of that
 *    boundary. (Same construction as papilo's Num::hashCode and HiGHS's
 *    double_hash_code, both of which document the 0.5 case at length.)
 *
 * The hash is a filter, never a proof: every candidate pair is verified
 * coefficient by coefficient before anything is deleted.
 *
 * ### Removes rows, so: MIP only
 *
 * A deleted row has no dual value, and no attempt is made to recover one --
 * there is no postsolve here. Cbc does not report duals for a MIP, so this is
 * free on the branch-and-bound path; it is *not* safe for the LP-only commands
 * (-dualSimplex, -primalSimplex, -barrier, -solveContinuous), which do report
 * them. Hence CbcSolver::preRootLPStrenghtening() runs this step only when its
 * caller passes allowRowRemoval=true, which only the BAB pipeline does.
 *
 * ### Usage
 * \code
 *   CbcRowReductions rr;
 *   if (!rr.run(solver, handler, logLevel))
 *     ; // either nothing to do, or infeasible -- check isInfeasible()
 *   printf("%d fixed, %d duplicate, %d parallel\n",
 *          rr.nFixedRows(), rr.nDuplicateRows(), rr.nParallelRows());
 * \endcode
 */
class CBCLIB_EXPORT CbcRowReductions {
public:
  CbcRowReductions();

  /*! \brief Remove redundant rows of \p solver in place.
   *
   * \param solver   The (already-loaded) MIP solver interface.
   * \param handler  Message handler (may be null); only used for logging.
   * \param logLevel CbcModel log level; >= 3 prints one line per removal.
   *
   * \return true if any row was removed or any row bound was tightened.
   *         Returns false both when there was nothing to do and when
   *         infeasibility was proved -- call isInfeasible() to distinguish.
   *         Nothing is deleted once infeasibility is known.
   */
  bool run(OsiSolverInterface *solver,
    CoinMessageHandler *handler,
    int logLevel);

  /// True if the reductions proved the model infeasible.
  bool isInfeasible() const { return infeasible_; }

  /// Rows removed because all their columns were fixed.
  int nFixedRows() const { return nFixedRows_; }

  /// Rows removed as exact duplicates of another row (lambda == 1).
  int nDuplicateRows() const { return nDuplicateRows_; }

  /// Rows removed as scalar multiples of another row (lambda != 1).
  int nParallelRows() const { return nParallelRows_; }

  /// Row bounds tightened while merging a duplicate/parallel pair.
  int nBoundsTightened() const { return nBoundsTightened_; }

  /// Total rows removed: nFixedRows() + nDuplicateRows() + nParallelRows().
  int nRowsRemoved() const
  {
    return nFixedRows_ + nDuplicateRows_ + nParallelRows_;
  }

  /// Wall-clock seconds used by run().
  double timeUsed() const { return timeUsed_; }

private:
  bool infeasible_;
  int nFixedRows_;
  int nDuplicateRows_;
  int nParallelRows_;
  int nBoundsTightened_;
  double timeUsed_;
};

#endif // CbcRowReductions_hpp
