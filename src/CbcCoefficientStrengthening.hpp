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

private:
  int nCoefficients_;
  int nRows_;
  double timeUsed_;
};

#endif // CbcCoefficientStrengthening_hpp
