// Copyright (C) 2026 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#ifndef CbcCoefficientStrengthening_hpp
#define CbcCoefficientStrengthening_hpp

#include "CbcConfig.h"

class OsiSolverInterface;
class CoinMessageHandler;

/*! \brief Coefficient tightening ("big-M strengthening") for the pre-root-LP
 *         strengthening phase.
 *
 * Shrinks oversized integer coefficients in one-sided rows, with a
 * compensating right-hand-side adjustment, so that the LP relaxation is
 * tighter while the *integer*-feasible set is preserved exactly. Like bound
 * propagation and clique merging, it belongs before the root LP solve because
 * it needs no LP information (no duals, no reduced costs) and removes no
 * column, so user callbacks still see their own variables.
 *
 * ### The rule
 *
 * Normalise a row with exactly one finite side to \f$a \cdot x \le b\f$ and
 * let \f$M\f$ be its maximum activity over the current column bounds. With
 * \f$s = M - b > 0\f$ (the row's slack; a non-positive slack means the row can
 * never be violated), any *integer* column \f$j\f$ with \f$|a_j| > s\f$ can be
 * shrunk to \f$\pm s\f$:
 *
 *  - \f$a_j > 0\f$: \f$a_j \leftarrow s\f$, \f$b \leftarrow b - (a_j - s)u_j\f$
 *  - \f$a_j < 0\f$: \f$a_j \leftarrow -s\f$, \f$b \leftarrow b - (a_j + s)l_j\f$
 *
 * At \f$x_j\f$ = its max-activity bound the two rows are algebraically
 * identical; one integer step away from it, *both* are redundant — which is
 * precisely where \f$|a_j| > s\f$ is needed, and why \f$x_j\f$ must be
 * integer. The slack \f$s\f$ is invariant under the change, so every
 * coefficient of a row uses the same \f$s\f$ and the order does not matter.
 *
 * A single pass is enough: because \f$s\f$ is invariant, a row cannot become
 * tightenable again, and tightening does not enable further bound propagation
 * (measured: identical fixpoint tightening counts before and after).
 *
 * ### Conflict-aware slack (clique cover)
 *
 * \f$M\f$ above assumes every unit-coefficient binary column of the row can
 * independently reach its upper bound at the same time. When the row also
 * has \f$k\ge 2\f$ such columns (coefficient exactly \f$+1\f$ after
 * normalisation) that pairwise conflict with each other (per the model's
 * conflict graph, see CoinConflictGraph -- typically populated by clique
 * merging earlier in the same pre-root-LP phase), that assumption is
 * pessimistic: at most one column of a mutual-conflict clique can be 1 at
 * once, so covering those \f$k\f$ columns with \f$p \le k\f$ cliques (a
 * greedy cover -- exact minimum clique cover is NP-hard) shows the true
 * maximum activity is \f$M - (k - p)\f$, not \f$M\f$. Substituting this
 * tighter bound for \f$M\f$ (equivalently, \f$s\f$ shrinks by \f$k-p\f$)
 * feeds directly into the same per-column rule above, so it can additionally
 * shrink *other* (non-clique) coefficients in the row -- most notably a
 * companion big-\f$M\f$ column, e.g. \f$\sum x_i \le M y\f$ reformulations,
 * where this often lowers \f$M\f$ well below the naive \f$k\f$.
 *
 * This is skipped (silently, not a correctness issue) whenever:
 *  - no conflict graph has been built yet (OsiSolverInterface::getCGraph()
 *    is null -- e.g. -clqStrengthening=off/after);
 *  - the row's own data alone would already force a pairwise conflict
 *    between two of those columns (guards against using a row to "improve"
 *    itself off conflicts that only exist because of the row itself --
 *    see the CLIQUECOVER_CIRCULARITY_TOL comment in the .cpp);
 *  - the row has more than \link CLIQUECOVER_MAX_ROW_VARS \endlink such
 *    columns, or the running total time spent on this sub-step across the
 *    whole \c run() call has passed \link CLIQUECOVER_TIME_BUDGET \endlink
 *    -- a naive-worst-case-quadratic greedy cover on a row with tens of
 *    thousands of columns can take seconds; both caps keep \c run() itself
 *    cheap on every corpus instance measured so far (see the class's git
 *    history for the benchmark).
 *
 * ### Usage
 * \code
 *   CbcCoefficientStrengthening cs;
 *   cs.run(solver, handler, logLevel);   // applied to solver in place
 *   printf("%d coefficients on %d rows\n", cs.nCoefficients(), cs.nRows());
 * \endcode
 */
class CBCLIB_EXPORT CbcCoefficientStrengthening {
public:
  CbcCoefficientStrengthening();

  /*! \brief Tighten coefficients of \p solver in place.
   *
   * \param solver   The (already-loaded) MIP solver interface.
   * \param handler  Message handler (may be null); only used for logging.
   * \param logLevel CbcModel log level; >= 3 prints one line per change.
   *
   * \return true if at least one coefficient was changed.
   *
   * Nothing is changed unless the solver actually honours
   * OsiSolverInterface::replaceMatrixOptional() -- the base class declares it
   * as a no-op, so for a solver that does not implement it this returns false
   * and leaves the model exactly as it was (row bounds included).
   */
  bool run(OsiSolverInterface *solver,
    CoinMessageHandler *handler,
    int logLevel);

  /// Number of coefficients whose magnitude was reduced.
  int nCoefficients() const { return nCoefficients_; }

  /// Number of rows in which at least one coefficient was reduced.
  int nRows() const { return nRows_; }

  /// Wall-clock seconds used by run().
  double timeUsed() const { return timeUsed_; }

  /// Number of rows where a clique cover found a tighter-than-naive slack
  /// (see the class comment's "Conflict-aware slack" section).
  int nCliqueCoverRows() const { return nCliqueCoverRows_; }

  /// Total slack reduction (sum of k - p over cliqueCoverRows) contributed
  /// by clique covers, i.e. how much tighter the naive per-row big-M bound
  /// became before the usual per-column shrink rule was applied.
  int cliqueCoverReduction() const { return cliqueCoverReduction_; }

  /// Wall-clock seconds spent computing clique covers (subset of timeUsed()).
  double cliqueCoverTimeUsed() const { return cliqueCoverTimeUsed_; }

private:
  int nCoefficients_;
  int nRows_;
  double timeUsed_;
  int nCliqueCoverRows_;
  int cliqueCoverReduction_;
  double cliqueCoverTimeUsed_;
};

#endif // CbcCoefficientStrengthening_hpp
