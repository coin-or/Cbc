// Copyright (C) 2024, MIPster contributors.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CoinPresolveParity_H
#define CoinPresolveParity_H

#include "CoinPresolveMatrix.hpp"

/*! \file

  GF(2) parity presolve action.

  Detects blocks of equality constraints that encode a binary parity
  (XOR / "Lights Out") subsystem and solves them exactly using bitwise
  Gaussian elimination over GF(2).

  A row qualifies as a *parity row* when:
    - rowLower == rowUpper  (equality),
    - every variable in the row is integer (typically binary),
    - all coefficients are in {1, -1, 2, -2},
    - at most one variable has |coefficient| == 2 (the "mod-2 helper" y_i).

  When a sufficient block of such rows is found (configurable threshold,
  default 10), the parity subsystem is extracted into a dense binary
  matrix packed in uint64_t words and solved via XOR-based row reduction.

  Possible outcomes:
    1. **Inconsistency** (0 = 1 in GF(2)) → declare the MIP infeasible.
    2. **Unique solution** → fix all participating binaries and helpers;
       remove the parity rows as redundant.
    3. **Free variables** → express dependent variables in terms of
       free ones, fix the determined binaries, and remove the
       consumed rows.
*/

/*! \class parity_action
    \brief Solve parity (XOR) subsystems via GF(2) Gaussian elimination.

  The presolve side detects parity rows, extracts the GF(2) system,
  performs bitwise Gaussian elimination, fixes determined variables,
  and removes consumed rows.

  The postsolve side restores the removed rows and variable bounds.
*/
class COINUTILSLIB_EXPORT parity_action : public CoinPresolveAction
{

public:
  /// Information needed to undo one fixed variable during postsolve
  struct fixed_var {
    int col;        ///< Column index of the fixed variable
    double origLo;  ///< Original lower bound
    double origUp;  ///< Original upper bound
    double fixVal;  ///< Value it was fixed to
  };

  /// Information needed to undo one removed row during postsolve
  struct removed_row {
    int row;          ///< Row index
    double rlo;       ///< Original row lower bound
    double rup;       ///< Original row upper bound
    int ninrow;       ///< Number of coefficients in the row
    const int *cols;  ///< Column indices (owned, allocated with new[])
    const double *els; ///< Coefficients (owned, allocated with new[])
  };

private:
  int nFixedVars_;                  ///< Number of fixed variables
  const fixed_var *fixedVars_;      ///< Array of fixed variable records
  int nRemovedRows_;                ///< Number of removed rows
  const removed_row *removedRows_;  ///< Array of removed row records

  /// Private constructor — only presolve() creates instances.
  parity_action(int nFixedVars,
    const fixed_var *fixedVars,
    int nRemovedRows,
    const removed_row *removedRows,
    const CoinPresolveAction *next)
    : CoinPresolveAction(next)
    , nFixedVars_(nFixedVars)
    , fixedVars_(fixedVars)
    , nRemovedRows_(nRemovedRows)
    , removedRows_(removedRows)
  {
  }

public:
  /// Returns "parity_action".
  const char *name() const;

  /*! \brief Detect and solve GF(2) parity subsystems.

    Scans the matrix for parity rows, extracts the GF(2) system,
    solves it, fixes variables, and removes consumed rows.

    \param prob  The presolve matrix to transform.
    \param next  Current head of the postsolve action list.
    \return New head of the postsolve action list (may be unchanged
            if no parity block was found).
  */
  static const CoinPresolveAction *presolve(CoinPresolveMatrix *prob,
    const CoinPresolveAction *next);

  /*! \brief Postsolve: restore removed rows and variable bounds. */
  void postsolve(CoinPostsolveMatrix *prob) const;

  /// Destructor
  virtual ~parity_action();
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
