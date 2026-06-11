// Copyright (C) 2024, MIPster contributors.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "CoinPresolveParity.hpp"
#include "CoinPresolveMatrix.hpp"
#include "CoinPresolveEmpty.hpp" // for DROP_COL/DROP_ROW
#include "CoinPresolveFixed.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"
#include "CoinMessage.hpp"

#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
#include "CoinPresolvePsdebug.hpp"
#endif

/*
  ================================================================
  GF(2) Dense Binary Matrix

  Rows of a GF(2) system are stored as arrays of uint64_t words.
  An augmented column (the RHS) is stored as a separate bit-vector.
  Arithmetic is XOR-based; division is trivial (1/1 = 1 in GF(2)).
  ================================================================
*/

namespace {

/// Number of uint64_t words needed for nCols bits
inline int wordsNeeded(int nCols)
{
  return (nCols + 63) / 64;
}

/// Get bit j in row
inline int getBit(const uint64_t *row, int j)
{
  return (row[j / 64] >> (j % 64)) & 1;
}

/// Set bit j in row to 1
inline void setBit(uint64_t *row, int j)
{
  row[j / 64] |= (uint64_t(1) << (j % 64));
}

/// XOR row src into row dst (in-place)
inline void xorRow(uint64_t *dst, const uint64_t *src, int nWords)
{
  for (int w = 0; w < nWords; w++)
    dst[w] ^= src[w];
}

/// Check if a row is entirely zero (excluding augmented column)
inline bool isZeroRow(const uint64_t *row, int nWords)
{
  for (int w = 0; w < nWords; w++)
    if (row[w] != 0)
      return false;
  return true;
}

/*
  ================================================================
  Phase 2: Fast Detection — O(nrows) Scan

  A row is a "parity candidate" if:
    (a) rlo[i] == rup[i]  (equality constraint)
    (b) every variable is integer
    (c) coefficients ∈ {1, -1, 2, -2}
    (d) at most one variable has |coefficient| == 2

  We also record, for each qualifying row, which column is the
  "helper" variable (|coeff| == 2) — or -1 if none.
  ================================================================
*/

struct ParityRowInfo {
  int row;        ///< Row index in the original matrix
  int helperCol;  ///< Column with |coeff| == 2, or -1
};

/// Scan matrix for parity-candidate rows.  Returns the list and count.
static void detectParityRows(
  CoinPresolveMatrix *prob,
  std::vector<ParityRowInfo> &candidates)
{
  const int nrows = prob->nrows_;
  const CoinBigIndex *mrstrt = prob->mrstrt_;
  const int *hinrow = prob->hinrow_;
  const double *rowels = prob->rowels_;
  const int *hcol = prob->hcol_;
  const double *rlo = prob->rlo_;
  const double *rup = prob->rup_;
  const double tol = 1.0e-8;

  candidates.clear();

  for (int i = 0; i < nrows; i++) {
    // Must be non-empty equality
    if (hinrow[i] == 0)
      continue;
    if (fabs(rlo[i] - rup[i]) > tol)
      continue;
    // Prohibited?
    if (prob->rowProhibited2(i))
      continue;

    const CoinBigIndex krs = mrstrt[i];
    const CoinBigIndex kre = krs + hinrow[i];

    bool valid = true;
    int helperCol = -1;
    int nHelpers = 0;

    for (CoinBigIndex k = krs; k < kre; k++) {
      int j = hcol[k];
      double aij = rowels[k];
      double absAij = fabs(aij);

      // All variables must be integer
      if (!prob->isInteger(j)) {
        valid = false;
        break;
      }
      // Coefficient must be 1, -1, 2, or -2
      if (fabs(absAij - 1.0) > tol && fabs(absAij - 2.0) > tol) {
        valid = false;
        break;
      }
      // Track helpers (|coeff| == 2)
      if (fabs(absAij - 2.0) <= tol) {
        nHelpers++;
        helperCol = j;
        if (nHelpers > 1) {
          valid = false;
          break;
        }
      }
      // Check variable bounds consistent with binary interpretation
      // (coefficient-1 variables should be binary for GF(2) to be exact)
      if (fabs(absAij - 1.0) <= tol) {
        // Variable should be binary (0/1 bounds)
        if (prob->clo_[j] < -tol || prob->cup_[j] > 1.0 + tol) {
          valid = false;
          break;
        }
      }
    }

    if (valid && nHelpers <= 1) {
      ParityRowInfo info;
      info.row = i;
      info.helperCol = helperCol;
      candidates.push_back(info);
      // Temporary debug print removed for cleanliness
    }
  }
}

/*
  ================================================================
  Phase 3: GF(2) Extraction

  Build a dense augmented binary matrix [A|b] over GF(2).

  For each qualifying row:
    - Variables with |coeff| == 1 contribute their column bit.
    - Variables with |coeff| == 2 are "helpers" (they don't appear
      in the GF(2) system because 2 ≡ 0 mod 2).
    - The RHS b_i = rlo[i] mod 2 (adjusted for helper contributions
      if there are fixed helpers, though typically they're free).

  We build a mapping:  gf2Col -> originalCol (for coeff-1 vars)
  ================================================================
*/

struct GF2System {
  int nRows;
  int nCols;          ///< Number of GF(2) columns (coeff-1 binary vars)
  int nWords;         ///< uint64_t words per row
  uint64_t *matrix;   ///< Dense GF(2) matrix, nRows × nWords
  uint64_t *rhs;      ///< Augmented column, nRows bits (one word per row, bit 0)
  int *gf2ColToOrig;  ///< gf2Col -> original column index
  int *origToGf2Col;  ///< original column -> gf2Col (or -1)
  std::vector<ParityRowInfo> rowInfo; ///< Per-row info (original row, helper col)

  GF2System()
    : nRows(0)
    , nCols(0)
    , nWords(0)
    , matrix(nullptr)
    , rhs(nullptr)
    , gf2ColToOrig(nullptr)
    , origToGf2Col(nullptr)
  {
  }

  ~GF2System()
  {
    delete[] matrix;
    delete[] rhs;
    delete[] gf2ColToOrig;
    delete[] origToGf2Col;
  }
};

/// Build the GF(2) system from the candidate rows.
static bool buildGF2System(
  CoinPresolveMatrix *prob,
  const std::vector<ParityRowInfo> &candidates,
  GF2System &sys)
{
  const CoinBigIndex *mrstrt = prob->mrstrt_;
  const int *hinrow = prob->hinrow_;
  const double *rowels = prob->rowels_;
  const int *hcol = prob->hcol_;
  const double *rlo = prob->rlo_;
  const double tol = 1.0e-8;

  // Step 1: Collect all binary (coeff-1) columns that appear
  //         in the candidate rows.
  int ncols = prob->ncols_;
  sys.origToGf2Col = new int[ncols];
  for (int j = 0; j < ncols; j++)
    sys.origToGf2Col[j] = -1;

  // Temporary vector to collect column indices
  std::vector<int> gf2Cols;
  for (int c = 0; c < static_cast<int>(candidates.size()); c++) {
    int i = candidates[c].row;
    CoinBigIndex krs = mrstrt[i];
    CoinBigIndex kre = krs + hinrow[i];
    for (CoinBigIndex k = krs; k < kre; k++) {
      int j = hcol[k];
      double absAij = fabs(rowels[k]);
      // Only coeff-1 variables go into GF(2)
      if (fabs(absAij - 1.0) <= tol && sys.origToGf2Col[j] == -1) {
        sys.origToGf2Col[j] = static_cast<int>(gf2Cols.size());
        gf2Cols.push_back(j);
      }
    }
  }

  sys.nCols = static_cast<int>(gf2Cols.size());
  sys.nRows = static_cast<int>(candidates.size());
  sys.nWords = wordsNeeded(sys.nCols);
  sys.rowInfo = candidates;

  if (sys.nCols == 0 || sys.nRows == 0)
    return false;

  // Build the gf2Col -> original mapping
  sys.gf2ColToOrig = new int[sys.nCols];
  for (int g = 0; g < sys.nCols; g++)
    sys.gf2ColToOrig[g] = gf2Cols[g];

  // Step 2: Allocate and fill the dense binary matrix
  int totalWords = sys.nRows * sys.nWords;
  sys.matrix = new uint64_t[totalWords];
  memset(sys.matrix, 0, totalWords * sizeof(uint64_t));
  sys.rhs = new uint64_t[sys.nRows];
  memset(sys.rhs, 0, sys.nRows * sizeof(uint64_t));

  for (int r = 0; r < sys.nRows; r++) {
    uint64_t *row = sys.matrix + r * sys.nWords;
    int i = candidates[r].row;
    CoinBigIndex krs = mrstrt[i];
    CoinBigIndex kre = krs + hinrow[i];

    // RHS = floor(rlo[i]) mod 2
    // (rlo == rup for equality; we take the integer value mod 2)
    int rhsBit = (static_cast<int>(fabs(rlo[i]) + 0.5)) % 2;
    // Handle negative RHS: -b ≡ b (mod 2) for odd b, 0 for even
    if (rlo[i] < -tol) {
      rhsBit = (static_cast<int>(fabs(rlo[i]) + 0.5)) % 2;
    }
    sys.rhs[r] = rhsBit;

    for (CoinBigIndex k = krs; k < kre; k++) {
      int j = hcol[k];
      double absAij = fabs(rowels[k]);
      if (fabs(absAij - 1.0) <= tol) {
        int gf2Col = sys.origToGf2Col[j];
        if (gf2Col >= 0)
          setBit(row, gf2Col);
      }
      // |coeff| == 2 variables are invisible in GF(2)
      // but if they have a contribution to the RHS we need to account
      // for it.  Since 2*y ≡ 0 (mod 2) for any integer y, the helper
      // variable does not affect the RHS parity.
    }
  }

  return true;
}

/*
  ================================================================
  Phase 4: Bitwise Gaussian Elimination over GF(2)

  Standard row reduction, but with XOR instead of floating-point ops.
  Each row elimination processes 64 columns per CPU instruction.

  Returns:
    0 = inconsistent (infeasible)
    1 = unique solution (all pivots found)
    2 = free variables exist (under-determined)

  On return, pivotCol[r] = GF(2) column index of the pivot in row r,
  or -1 if row r has no pivot (either zero row or inconsistent).
  solution[g] = 0 or 1 for each GF(2) column g.
  ================================================================
*/

static int gaussianEliminationGF2(
  GF2System &sys,
  int *pivotCol,
  int *solution)
{
  int nRows = sys.nRows;
  int nCols = sys.nCols;
  int nWords = sys.nWords;
  uint64_t *matrix = sys.matrix;
  uint64_t *rhs = sys.rhs;

  // Initialize pivot tracking
  for (int r = 0; r < nRows; r++)
    pivotCol[r] = -1;

  // Forward elimination (row echelon form)
  int pivotRow = 0;
  int rank = 0;

  for (int col = 0; col < nCols && pivotRow < nRows; col++) {
    // Find a row with a 1 in this column, at or below pivotRow
    int found = -1;
    for (int r = pivotRow; r < nRows; r++) {
      uint64_t *row = matrix + r * nWords;
      if (getBit(row, col)) {
        found = r;
        break;
      }
    }
    if (found < 0)
      continue; // No pivot in this column; it's a free variable

    // Swap found row with pivotRow
    if (found != pivotRow) {
      uint64_t *rowA = matrix + pivotRow * nWords;
      uint64_t *rowB = matrix + found * nWords;
      for (int w = 0; w < nWords; w++) {
        uint64_t tmp = rowA[w];
        rowA[w] = rowB[w];
        rowB[w] = tmp;
      }
      uint64_t tmpRhs = rhs[pivotRow];
      rhs[pivotRow] = rhs[found];
      rhs[found] = tmpRhs;
      // Also swap rowInfo so we track which original row is where
      std::swap(sys.rowInfo[pivotRow], sys.rowInfo[found]);
    }

    pivotCol[pivotRow] = col;

    // Eliminate this column from all other rows
    uint64_t *pivRow = matrix + pivotRow * nWords;
    for (int r = 0; r < nRows; r++) {
      if (r == pivotRow)
        continue;
      uint64_t *row = matrix + r * nWords;
      if (getBit(row, col)) {
        xorRow(row, pivRow, nWords);
        rhs[r] ^= rhs[pivotRow];
      }
    }

    pivotRow++;
    rank++;
  }

  // Check for inconsistency: any row that is all-zero in the
  // matrix but has rhs == 1 means 0 = 1, i.e., infeasible.
  for (int r = 0; r < nRows; r++) {
    uint64_t *row = matrix + r * nWords;
    if (isZeroRow(row, nWords) && rhs[r] != 0) {
#if PRESOLVE_DEBUG > 0
      printf("PARITY_DEBUG: GF(2) Inconsistency found at pivot-row %d (original row %d), rhs=%d\n",
             r, sys.rowInfo[r].row, (int)rhs[r]);
#endif
      return 0; // Infeasible
    }
  }

  // Back-substitution to find solution
  // Initialize all solutions to 0 (free variables default to 0)
  for (int g = 0; g < nCols; g++)
    solution[g] = 0;

  // Process pivot rows in reverse order
  for (int r = rank - 1; r >= 0; r--) {
    int pc = pivotCol[r];
    if (pc < 0)
      continue;
    // solution[pc] = rhs[r] XOR sum of (matrix[r][g] * solution[g]) for g != pc
    int val = static_cast<int>(rhs[r]);
    uint64_t *row = matrix + r * nWords;
    for (int g = 0; g < nCols; g++) {
      if (g != pc && getBit(row, g))
        val ^= solution[g];
    }
    solution[pc] = val;
  }

  return (rank == nCols) ? 1 : 2;
}

} /* end unnamed namespace */

/*
  ================================================================
  Phase 5 & 6: The Presolve Action
  ================================================================
*/

const char *parity_action::name() const
{
  return "parity_action";
}

const CoinPresolveAction *parity_action::presolve(
  CoinPresolveMatrix *prob,
  const CoinPresolveAction *next)
{
#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
#if PRESOLVE_DEBUG > 0
  std::cout
    << "Entering parity_action::presolve, "
    << prob->nrows_ << " rows, "
    << prob->ncols_ << " cols." << std::endl;
#endif
  presolve_check_sol(prob);
  presolve_check_nbasic(prob);
#endif

  // Phase 2: Detect parity rows
  std::vector<ParityRowInfo> candidates;
  detectParityRows(prob, candidates);

  // Bail out if fewer than 10 parity rows found
  const int minParityRows = 10;
  if (static_cast<int>(candidates.size()) < minParityRows) {
#if PRESOLVE_DEBUG > 0
    std::cout
      << "parity_action: only " << candidates.size()
      << " parity rows found (need " << minParityRows
      << "), skipping." << std::endl;
#endif
    return next;
  }

  // Phase 3: Build GF(2) system
  GF2System sys;
  if (!buildGF2System(prob, candidates, sys)) {
#if PRESOLVE_DEBUG > 0
    printf("PARITY_DEBUG: failed to build GF(2) system.\n");
#endif
    return next;
  }

#if PRESOLVE_DEBUG > 0
  printf("PARITY_DEBUG: GF(2) system has %d rows, %d cols.\n", sys.nRows, sys.nCols);
#endif

  // Phase 4: Gaussian elimination
  int *pivotCol = new int[sys.nRows];
  int *solution = new int[sys.nCols];
  int result = gaussianEliminationGF2(sys, pivotCol, solution);

  if (result == 0) {
    // Infeasible!
    prob->status_ = 1;
    delete[] pivotCol;
    delete[] solution;
#if PRESOLVE_DEBUG > 0
    printf("PARITY_DEBUG: GF(2) system is inconsistent -> INFEASIBLE.\n");
#endif
    return next;
  }

#if PRESOLVE_DEBUG > 0
  printf("PARITY_DEBUG: GF(2) elimination result = %s, rank = ",
         (result == 1 ? "unique" : "free variables"));
  {
    int rank = 0;
    for (int r = 0; r < sys.nRows; r++)
      if (pivotCol[r] >= 0)
        rank++;
    printf("%d", rank);
  }
  printf(" / %d\n", sys.nCols);
#endif

  // Phase 5: Fix variables and remove rows
  //
  // We fix *all* determined binary variables (those that are pivot columns).
  // Free variables (non-pivot columns) are left alone with solution[g]=0,
  // but we don't fix them since other values may also satisfy the parity.
  //
  // For the "helper" variables (|coeff| == 2), we back-substitute:
  //   2*y_i + sum(±x_j) = b_i  →  y_i = (b_i - sum(±x_j*solution_j)) / 2

  // Collect variable fixings
  std::vector<fixed_var> fixings;
  // Collect row removals
  std::vector<removed_row> removals;

  // Step 5a: Fix determined binaries (pivot columns only)
  // We only fix variables whose GF(2) column was a pivot.
  // If the system is under-determined (result == 2), we can only fix a pivot
  // if it is independent of all free variables (i.e. its row has exactly one 1).
  for (int r = 0; r < sys.nRows; r++) {
    int pc = pivotCol[r];
    if (pc < 0)
      continue;
      
    if (result == 2) {
      int bitCount = 0;
      uint64_t *row = sys.matrix + r * sys.nWords;
      for (int w = 0; w < sys.nWords; w++) {
        uint64_t word = row[w];
        while (word) {
          bitCount++;
          word &= (word - 1);
        }
      }
      if (bitCount > 1)
        continue; // Depends on free variables, cannot fix
    }

    int origCol = sys.gf2ColToOrig[pc];
    double val = static_cast<double>(solution[pc]);

    // Only fix if bounds permit (they should, since we checked binary)
    if (val >= prob->clo_[origCol] - 1.0e-8
      && val <= prob->cup_[origCol] + 1.0e-8
      && fabs(prob->clo_[origCol] - prob->cup_[origCol]) > 1.0e-8) {
      // Don't fix if prohibited
      if (prob->colProhibited2(origCol))
        continue;

      fixed_var fv;
      fv.col = origCol;
      fv.origLo = prob->clo_[origCol];
      fv.origUp = prob->cup_[origCol];
      fv.fixVal = val;
      fixings.push_back(fv);

      // Fix the variable
      prob->clo_[origCol] = val;
      prob->cup_[origCol] = val;
      if (prob->sol_)
        prob->sol_[origCol] = val;
      prob->addCol(origCol);
    }
  }

  if (fixings.empty()) {
    // Nothing to do — maybe all were prohibited
    delete[] pivotCol;
    delete[] solution;
    return next;
  }

  // Step 5b: For each parity row with a pivot, compute helper value
  //          and remove the row.
  const CoinBigIndex *mrstrt = prob->mrstrt_;
  int *hinrow = prob->hinrow_;
  const double *rowels = prob->rowels_;
  const int *hcol = prob->hcol_;
  double *rlo = prob->rlo_;
  double *rup = prob->rup_;
  double *colels = prob->colels_;
  int *hrow = prob->hrow_;
  const CoinBigIndex *mcstrt = prob->mcstrt_;
  int *hincol = prob->hincol_;
  const double tol = 1.0e-8;

  // We CANNOT do this if there are free variables (result == 2), because the
  // original rows still enforce constraints among the unfixed free variables.
  if (result == 1) {
    for (int r = 0; r < sys.nRows; r++) {
      // Only remove rows where we have a pivot (i.e., the row was
      // used in the elimination and is now redundant given the fixings)
      if (pivotCol[r] < 0)
        continue;

      int origRow = sys.rowInfo[r].row;
      int helperCol = sys.rowInfo[r].helperCol;

    // Don't remove prohibited rows
    if (prob->rowProhibited2(origRow))
      continue;

    // Compute helper variable value if present
    if (helperCol >= 0
      && !prob->colProhibited2(helperCol)
      && fabs(prob->clo_[helperCol] - prob->cup_[helperCol]) > tol) {
      // 2*y + sum(a_j * x_j) = b
      // For coeff-1 vars: contribution = sum(±1 * solution_j)
      // y = (b - contribution) / 2
      double b = rlo[origRow];
      double contribution = 0.0;
      CoinBigIndex krs = mrstrt[origRow];
      CoinBigIndex kre = krs + hinrow[origRow];
      double helperCoeff = 0.0;

      for (CoinBigIndex k = krs; k < kre; k++) {
        int j = hcol[k];
        double aij = rowels[k];
        if (j == helperCol) {
          helperCoeff = aij;
        } else {
          // This variable should now be fixed
          double xj = prob->clo_[j]; // fixed value
          contribution += aij * xj;
        }
      }

      if (fabs(helperCoeff) > tol) {
        double helperVal = (b - contribution) / helperCoeff;
        // Round to nearest integer (should be exact)
        helperVal = floor(helperVal + 0.5);

        // Fix the helper
        if (helperVal >= prob->clo_[helperCol] - tol
          && helperVal <= prob->cup_[helperCol] + tol) {
          fixed_var fv;
          fv.col = helperCol;
          fv.origLo = prob->clo_[helperCol];
          fv.origUp = prob->cup_[helperCol];
          fv.fixVal = helperVal;
          fixings.push_back(fv);

          prob->clo_[helperCol] = helperVal;
          prob->cup_[helperCol] = helperVal;
          if (prob->sol_)
            prob->sol_[helperCol] = helperVal;
          prob->addCol(helperCol);
        }
      }
    }

    // Save row for postsolve
    CoinBigIndex krs = mrstrt[origRow];
    int len = hinrow[origRow];

    removed_row rr;
    rr.row = origRow;
    rr.rlo = rlo[origRow];
    rr.rup = rup[origRow];
    rr.ninrow = len;
    rr.cols = CoinCopyOfArray(&hcol[krs], len);
    rr.els = CoinCopyOfArray(&rowels[krs], len);
    removals.push_back(rr);

    // Remove row from matrix
    // Delete this row from each column's entry list
    CoinBigIndex kre = krs + len;
    for (CoinBigIndex k = krs; k < kre; k++) {
      presolve_delete_from_col(origRow, hcol[k], mcstrt, hincol, hrow, colels);
      if (hincol[hcol[k]] == 0) {
        PRESOLVE_REMOVE_LINK(prob->clink_, hcol[k]);
      }
    }
    hinrow[origRow] = 0;
    PRESOLVE_REMOVE_LINK(prob->rlink_, origRow);

      // Zero out bounds (standard for removed rows)
      rlo[origRow] = 0.0;
      rup[origRow] = 0.0;
    }
  }

  delete[] pivotCol;
  delete[] solution;

  if (removals.empty() && fixings.empty()) {
    return next;
  }

  // Allocate and copy the action data
  int nFix = static_cast<int>(fixings.size());
  int nRem = static_cast<int>(removals.size());

  fixed_var *fvArr = nullptr;
  if (nFix > 0) {
    fvArr = new fixed_var[nFix];
    for (int i = 0; i < nFix; i++)
      fvArr[i] = fixings[i];
  }

  removed_row *rrArr = nullptr;
  if (nRem > 0) {
    rrArr = new removed_row[nRem];
    for (int i = 0; i < nRem; i++)
      rrArr[i] = removals[i];
  }

  CoinMessages messages = CoinMessage(prob->messages().language());
  prob->messageHandler()->message(COIN_PRESOLVE_PASS, messages)
    << nFix << nRem << -1 << CoinMessageEol;

#if PRESOLVE_DEBUG > 0
  std::cout
    << "parity_action: fixed " << nFix
    << " variables, removed " << nRem
    << " rows." << std::endl;
#endif

#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  presolve_check_sol(prob);
  presolve_check_nbasic(prob);
#endif

  return new parity_action(nFix, fvArr, nRem, rrArr, next);
}

/*
  ================================================================
  Phase 6: Postsolve — Reconstruct the original problem
  ================================================================
*/

void parity_action::postsolve(CoinPostsolveMatrix *prob) const
{
#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
#if PRESOLVE_DEBUG > 0
  std::cout
    << "Entering parity_action::postsolve, "
    << nFixedVars_ << " vars, "
    << nRemovedRows_ << " rows." << std::endl;
#endif
  presolve_check_sol(prob);
  presolve_check_nbasic(prob);
#endif

  double *colels = prob->colels_;
  int *hrow = prob->hrow_;
  CoinBigIndex *mcstrt = prob->mcstrt_;
  CoinBigIndex *link = prob->link_;
  int *hincol = prob->hincol_;

  double *rlo = prob->rlo_;
  double *rup = prob->rup_;
  double *sol = prob->sol_;
  double *acts = prob->acts_;
  double *clo = prob->clo_;
  double *cup = prob->cup_;

  CoinBigIndex &free_list = prob->free_list_;

  // Step 1: Restore removed rows (in reverse order)
  for (int i = nRemovedRows_ - 1; i >= 0; i--) {
    const removed_row &rr = removedRows_[i];

    rlo[rr.row] = rr.rlo;
    rup[rr.row] = rr.rup;

    double rowact = 0.0;
    for (int k = 0; k < rr.ninrow; k++) {
      int jcol = rr.cols[k];

      // Re-insert coefficient into threaded column storage
      {
        CoinBigIndex kk = free_list;
        assert(kk >= 0 && kk < prob->bulk0_);
        free_list = link[free_list];
        hrow[kk] = rr.row;
        colels[kk] = rr.els[k];
        link[kk] = mcstrt[jcol];
        mcstrt[jcol] = kk;
      }
      hincol[jcol]++;

      rowact += rr.els[k] * sol[jcol];
    }
#if PRESOLVE_CONSISTENCY
    presolve_check_free_list(prob);
#endif

    acts[rr.row] = rowact;

    // The row was basic (slack was basic) before removal; restore that.
    prob->setRowStatus(rr.row, CoinPrePostsolveMatrix::basic);
  }

  // Step 2: Restore variable bounds (in reverse order)
  for (int i = nFixedVars_ - 1; i >= 0; i--) {
    const fixed_var &fv = fixedVars_[i];
    clo[fv.col] = fv.origLo;
    cup[fv.col] = fv.origUp;
    // Variable value stays at fv.fixVal (it was correctly fixed)
    // Set status to reflect the fixed value vs bounds
    if (fabs(sol[fv.col] - fv.origLo) <= 1.0e-8) {
      prob->setColumnStatus(fv.col, CoinPrePostsolveMatrix::atLowerBound);
    } else if (fabs(sol[fv.col] - fv.origUp) <= 1.0e-8) {
      prob->setColumnStatus(fv.col, CoinPrePostsolveMatrix::atUpperBound);
    } else {
      // Value is between bounds — make it basic; we need to find
      // a row to make nonbasic in exchange. The first restored
      // parity row should serve.
      prob->setColumnStatus(fv.col, CoinPrePostsolveMatrix::basic);
    }
  }

#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  presolve_check_threads(prob);
  presolve_check_sol(prob);
  presolve_check_nbasic(prob);
#if PRESOLVE_DEBUG > 0
  std::cout << "Leaving parity_action::postsolve." << std::endl;
#endif
#endif
}

/*
  ================================================================
  Destructor
  ================================================================
*/

parity_action::~parity_action()
{
  for (int i = 0; i < nRemovedRows_; i++) {
    deleteAction(removedRows_[i].cols, int *);
    deleteAction(removedRows_[i].els, double *);
  }
  deleteAction(removedRows_, removed_row *);
  deleteAction(fixedVars_, fixed_var *);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
