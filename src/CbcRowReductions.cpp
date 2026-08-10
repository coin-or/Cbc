// Copyright (C) 2026 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcRowReductions.hpp"

#include "CoinMessageHandler.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

// Reference solution loaded by -debugCuts; set before applyLpMethod() so that
// pre-root-LP strengthening can check itself before the OsiRowCutDebugger is
// active. Same globals CbcBoundPropagation.cpp and
// CbcCoefficientStrengthening.cpp fall back on -- see the comment there for
// why the debugger itself is not yet attached this early when driven from the
// "-debugCuts" CLI switch.
extern double *debugSolution;
extern int debugNumberColumns;

namespace {

/// Anything at or beyond this magnitude is treated as infinite, matching the
/// convention used throughout Osi/Cgl (see e.g. CglPreProcess) and the sibling
/// CbcCoefficientStrengthening pass.
const double ROWRED_INFINITY = 1.0e20;

/// Feasibility tolerance, applied to row bounds and row activities. Same role
/// as papilo's Num::feastol and HiGHS's primal_feastol.
const double ROWRED_FEASTOL = 1.0e-6;

/// Bound-comparison epsilon: two bounds this close are the same bound.
const double ROWRED_EPSILON = 1.0e-9;

/// Coefficient-matching tolerance, applied *absolutely* to the difference
/// between a coefficient and its scaled counterpart. Deliberately much tighter
/// than ROWRED_FEASTOL: both papilo (Num::epsilon, 1e-9) and HiGHS
/// (small_matrix_value, 1e-9) hold coefficients to this while holding bounds
/// only to the feasibility tolerance, and conflating the two would let rows
/// that are merely *similar* be treated as parallel.
const double ROWRED_COEFTOL = 1.0e-9;

/// Refuse to merge a pair whose scale factor is this extreme. Such a pair is
/// mathematically fine but the transferred bound is then dominated by rounding
/// in the ratio, and losing the reduction costs nothing.
const double ROWRED_MAX_RATIO = 1.0e12;

/// How many surviving rows of a hash bucket a new row is compared against.
/// The O(k^2) blow-up on a bucket of k rows only happens when the bucket holds
/// many *distinct* parallel classes -- a single class of k mutually parallel
/// rows is fully collapsed by k-1 comparisons against one representative. So a
/// small cap bounds the worst case while still catching the common cases.
/// (papilo effectively uses 1 -- it compares only against bucket[0] -- and
/// notes the missing cap as a TODO; HiGHS has no cap at all.)
const int ROWRED_MAX_REPRESENTATIVES = 8;

/// splitmix64 finalizer: cheap, and mixes low bits into high ones so that
/// adding per-entry hashes together does not pile everything into one bucket.
inline uint64_t mix64(uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

/// Hash a *scale-normalised* coefficient into a coarse bucket, so that two
/// coefficients equal to within rounding noise land in the same bucket.
///
/// Only the exponent and the leading ~16 mantissa bits are kept -- more bits
/// would make the buckets so small that values equal within epsilon start
/// missing each other. The multiplication by 1/phi is what makes that work in
/// practice: without it a "round" value such as 0.5 sits exactly on a
/// power-of-two boundary, so 0.5 and 0.5-1e-16 get different exponents and
/// different buckets. Scaling by an irrational number moves the input off
/// those boundaries, and the values that land on one instead are irrational
/// multiples of a power of two, which do not occur in real model data. Both
/// papilo (Num::hashCode) and HiGHS (double_hash_code) use this same trick and
/// document the 0.5 case explicitly.
inline uint64_t bucketCoefficient(double value)
{
  if (value == 0.0)
    return 0;
  const double inverseGoldenRatio = 2.0 / (1.0 + std::sqrt(5.0));
  int exponent = 0;
  const double mantissa = std::frexp(value * inverseGoldenRatio, &exponent);
  // |mantissa| is in [0.5, 1), so this keeps 16 significant bits and the sign.
  const int64_t bits = static_cast< int64_t >(mantissa * 65536.0);
  return (static_cast< uint64_t >(static_cast< uint32_t >(bits)) << 32)
    ^ static_cast< uint64_t >(static_cast< uint32_t >(exponent));
}

/// Relative-difference test, papilo's Num::isSafeGT. Used *in addition to* the
/// absolute feasibility tolerance before declaring infeasibility, never for
/// declaring a row redundant: a wrong redundancy claim drops one constraint,
/// a wrong infeasibility claim throws away the whole model.
inline bool safeGreater(double a, double b)
{
  const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
  return (a - b) / scale > 1024.0 * std::numeric_limits< double >::epsilon();
}

/// One row, keyed for bucketing. Sorted by (hash, row) -- the row index as
/// final key keeps the outcome independent of anything but the model itself.
struct RowKey {
  uint64_t hash;
  int row;
};

inline bool operator<(const RowKey &a, const RowKey &b)
{
  if (a.hash != b.hash)
    return a.hash < b.hash;
  return a.row < b.row;
}

} // end anonymous namespace

CbcRowReductions::CbcRowReductions()
  : infeasible_(false)
  , nFixedRows_(0)
  , nDuplicateRows_(0)
  , nParallelRows_(0)
  , nBoundsTightened_(0)
  , timeUsed_(0.0)
{
}

bool CbcRowReductions::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel)
{
  infeasible_ = false;
  nFixedRows_ = 0;
  nDuplicateRows_ = 0;
  nParallelRows_ = 0;
  nBoundsTightened_ = 0;
  timeUsed_ = 0.0;

  if (!solver || solver->getNumRows() == 0)
    return false;

  const double t0 = CoinWallclockTime();

  const int numberRows = solver->getNumRows();
  const int numberColumns = solver->getNumCols();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();

  // Read-only throughout: every change is staged in the working vectors below
  // and applied to the solver in one batch at the very end. That keeps these
  // pointers valid for the whole scan (setRowLower()/deleteRows() would
  // invalidate the cached row copy), and means an infeasibility discovered
  // half-way through leaves the model exactly as it was.
  const CoinPackedMatrix *rowCopy = solver->getMatrixByRow();
  const double *element = rowCopy->getElements();
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();

  std::vector< char > deleted(numberRows, 0);
  std::vector< double > workLower(rowLower, rowLower + numberRows);
  std::vector< double > workUpper(rowUpper, rowUpper + numberRows);
  std::vector< char > boundChanged(numberRows, 0);

  // --- Step A: rows all of whose columns are fixed ------------------------
  //
  // Such a row has a constant activity, so it is either violated (the model is
  // infeasible) or it can never be binding (it goes). Rows that are only
  // *partly* fixed are deliberately left alone: folding a fixed column's
  // contribution into the row bounds and shrinking the row is bound
  // propagation's and CglPreProcess' job, and doing it here would rewrite row
  // content that CglPreProcess re-derives from scratch anyway.
  for (int iRow = 0; iRow < numberRows; iRow++) {
    const CoinBigIndex start = rowStart[iRow];
    const int length = rowLength[iRow];

    // Activity at the fixed values, plus a bound on how far it can move given
    // that the columns are only fixed to within a tolerance. This is the
    // per-row form of HiGHS' coefficient-aware fixed test: a column whose
    // bounds are 1e-9 apart is *not* effectively fixed if its coefficient is
    // 1e9, and testing the product per row is both cheaper and sharper than
    // testing max|a_j| per column.
    double activity = 0.0;
    double wobble = 0.0;
    bool allFixed = true;
    for (CoinBigIndex j = start; j < start + length; j++) {
      const int iColumn = column[j];
      const double spread = colUpper[iColumn] - colLower[iColumn];
      if (spread > ROWRED_EPSILON) {
        allFixed = false;
        break;
      }
      activity += element[j] * colLower[iColumn];
      wobble += std::fabs(element[j]) * spread;
    }
    // The wobble test here is an early-out, not a safety net: the satisfaction
    // test below already weighs the same wobble, so dropping this clause changes
    // no outcome (confirmed by mutation -- it is the coefficient *weighting*
    // that is load-bearing, and neutralising that is caught). It is kept because
    // it skips the bound loads and branches for a row that cannot qualify.
    if (!allFixed || wobble > ROWRED_FEASTOL)
      continue;

    const double lo = workLower[iRow];
    const double hi = workUpper[iRow];
    const bool hasLower = lo > -ROWRED_INFINITY;
    const bool hasUpper = hi < ROWRED_INFINITY;

    // Clearly violated: infeasible. Both an absolute and a relative test must
    // agree (see safeGreater) before we make that call.
    if ((hasUpper && (activity - wobble) - hi > ROWRED_FEASTOL
          && safeGreater(activity - wobble, hi))
      || (hasLower && lo - (activity + wobble) > ROWRED_FEASTOL
        && safeGreater(lo, activity + wobble))) {
      if (logLevel >= 1)
        printf("  Row reductions: row %d has all columns fixed and activity "
               "%.12g outside [%.12g,%.12g] — infeasible.\n",
          iRow, activity, lo, hi);
      infeasible_ = true;
      timeUsed_ = CoinWallclockTime() - t0;
      return false;
    }

    // Comfortably satisfied: the row can never bind, so drop it. Note the
    // threshold here is the *tight* tolerance, not the feasibility tolerance
    // used just above, and the gap between the two is deliberate: it leaves a
    // sliver -- a violation between 1e-9 and 1e-6 -- where the row is neither
    // declared infeasible nor deleted, but simply kept, so that the LP decides.
    //
    // Using the feasibility tolerance on both tests would close that sliver and
    // make every row that is not proven infeasible removable, which errs in the
    // one direction that matters: the row is *violated*, so deleting it turns an
    // infeasible model into a feasible-looking one, and it does so by up to
    // 1e-6, ten times Clp's own primal tolerance. Below 1e-9 the violation is
    // under the LP's noise floor and deleting the row changes nothing it could
    // have observed; a non-zero threshold rather than an exact test is what
    // keeps a satisfied equation from being kept over 1e-16 of rounding.
    const bool safeUpper = !hasUpper || (activity + wobble) - hi <= ROWRED_EPSILON;
    const bool safeLower = !hasLower || lo - (activity - wobble) <= ROWRED_EPSILON;
    if (safeUpper && safeLower) {
      deleted[iRow] = 1;
      nFixedRows_++;
      if (logLevel >= 3)
        printf("  Row reductions: row %d all columns fixed, activity %g in "
               "[%g,%g] — removed.\n",
          iRow, activity, lo, hi);
    }
  }

  // --- Step B: duplicate and parallel rows --------------------------------
  //
  // One mechanism for both: an exact duplicate is just a parallel row with
  // scale factor 1. Each row is normalised by its pivot before hashing, so
  // lambda*row and row hash identically and land in the same bucket; only rows
  // inside a bucket are ever compared, and every candidate pair is then
  // verified coefficient by coefficient. The hash is a filter, never a proof.
  std::vector< RowKey > keys;
  std::vector< double > pivotValue(numberRows, 0.0);
  std::vector< int > pivotColumn(numberRows, -1);
  keys.reserve(numberRows);

  for (int iRow = 0; iRow < numberRows; iRow++) {
    if (deleted[iRow])
      continue;
    const CoinBigIndex start = rowStart[iRow];
    const int length = rowLength[iRow];
    // A singleton row is a bound in disguise; bound propagation has already
    // run and owns that case. Both reference implementations skip these too.
    if (length <= 1)
      continue;

    // The pivot: largest magnitude, ties broken by *smallest column index*.
    // The tie-break is a correctness requirement, not a nicety -- two rows
    // related by lambda = -1 have equal-magnitude candidate pivots, and if
    // they picked different ones they would normalise to different vectors and
    // never meet in a bucket.
    // Comparison is deliberately *exact*, not within a tolerance, and this is
    // the more robust choice rather than a sloppier one: exact scaling of a row
    // preserves exact equality of magnitudes, so two rows related by
    // lambda = -1 still agree on which entries are tied and pick the same
    // column. A tolerance band would not scale with the row and could
    // therefore make the two disagree. Same comparison HiGHS uses.
    //
    // Note the tie-break cannot be exercised through OsiClp: Clp normalises
    // every row to ascending column order, so "first strict maximum in storage
    // order" already *is* "smallest column index among the tied maxima", and
    // mutating the tie-break away changes nothing measurable here. It is kept
    // because CoinPackedMatrix promises no index order and this pass takes an
    // OsiSolverInterface, so the guarantee is the interface's, not Clp's.
    int bestColumn = -1;
    double bestValue = 0.0;
    double bestMagnitude = -1.0;
    for (CoinBigIndex j = start; j < start + length; j++) {
      const double magnitude = std::fabs(element[j]);
      if (magnitude > bestMagnitude
        || (magnitude == bestMagnitude && column[j] < bestColumn)) {
        bestMagnitude = magnitude;
        bestColumn = column[j];
        bestValue = element[j];
      }
    }
    if (bestColumn < 0 || bestMagnitude <= ROWRED_COEFTOL)
      continue; // all-zero row: nothing to normalise by

    pivotValue[iRow] = bestValue;
    pivotColumn[iRow] = bestColumn;

    // Order-independent combination: the per-entry hashes are simply summed,
    // so the result does not depend on the order of the indices within the
    // row. CoinPackedMatrix makes no promise about that order, and this avoids
    // needing a sorted copy of every row (papilo can memcmp its index arrays
    // only because its own storage is sorted).
    uint64_t hash = mix64(static_cast< uint64_t >(length));
    for (CoinBigIndex j = start; j < start + length; j++) {
      const uint64_t entry = mix64(static_cast< uint64_t >(column[j]) * 2654435761ULL)
        ^ bucketCoefficient(element[j] / bestValue);
      hash += mix64(entry);
    }

    RowKey key;
    key.hash = hash;
    key.row = iRow;
    keys.push_back(key);
  }

  std::sort(keys.begin(), keys.end());

  // Scatter space for pair verification, allocated only once a bucket with
  // more than one row actually turns up. `mark` holds a generation stamp
  // rather than a boolean so it never needs clearing.
  std::vector< int > mark;
  std::vector< double > scatter;
  int generation = 0;

  const size_t numberKeys = keys.size();
  for (size_t bucketStart = 0; bucketStart < numberKeys;) {
    size_t bucketEnd = bucketStart + 1;
    while (bucketEnd < numberKeys && keys[bucketEnd].hash == keys[bucketStart].hash)
      bucketEnd++;
    if (bucketEnd - bucketStart < 2) {
      bucketStart = bucketEnd;
      continue;
    }

    if (mark.empty()) {
      mark.assign(numberColumns, -1);
      scatter.assign(numberColumns, 0.0);
    }

    // Rows of this bucket that have survived so far, in increasing row order.
    int representative[ROWRED_MAX_REPRESENTATIVES];
    int numberRepresentatives = 0;

    for (size_t k = bucketStart; k < bucketEnd; k++) {
      const int iRow = keys[k].row;
      if (deleted[iRow])
        continue;

      const CoinBigIndex startNew = rowStart[iRow];
      const int lengthNew = rowLength[iRow];

      // Scatter the new row once, then walk each representative against it.
      const int stamp = ++generation;
      for (CoinBigIndex j = startNew; j < startNew + lengthNew; j++) {
        mark[column[j]] = stamp;
        scatter[column[j]] = element[j];
      }

      bool merged = false;
      for (int rep = 0; rep < numberRepresentatives && !merged; rep++) {
        const int jRow = representative[rep];
        if (deleted[jRow])
          continue;
        const CoinBigIndex startRep = rowStart[jRow];
        const int lengthRep = rowLength[jRow];
        // Equal length plus "every column of the representative is in the new
        // row" gives equal support as a set. Both are re-checks rather than
        // discoveries -- the hash is seeded with the length, so unequal lengths
        // land in different buckets and are never compared, and a differing
        // support would have to forge a full hash collision to get here. They
        // cost two comparisons and mean a collision cannot delete a row.
        if (lengthRep != lengthNew)
          continue;

        // Scale factor between the two rows, read off the representative's
        // pivot column. Note this is the true ratio between the rows whether
        // or not the two pivots happen to sit in the same column.
        const int cPivot = pivotColumn[jRow];
        if (mark[cPivot] != stamp)
          continue;
        const double ratioNewOverRep = scatter[cPivot] / pivotValue[jRow];
        if (ratioNewOverRep == 0.0 || !std::isfinite(ratioNewOverRep))
          continue;
        if (std::fabs(ratioNewOverRep) > ROWRED_MAX_RATIO
          || std::fabs(ratioNewOverRep) < 1.0 / ROWRED_MAX_RATIO)
          continue;

        // Verify every coefficient, always scaling the *smaller* row up to the
        // larger one so that the absolute tolerance is applied to the larger
        // magnitudes. That direction can only reject a genuine pair, never
        // accept a false one -- the safe way round, since a false accept
        // deletes a real constraint.
        const bool newIsLarger = std::fabs(ratioNewOverRep) >= 1.0;
        const double up = newIsLarger ? ratioNewOverRep : 1.0 / ratioNewOverRep;
        bool parallel = true;
        for (CoinBigIndex j = startRep; j < startRep + lengthRep; j++) {
          const int iColumn = column[j];
          if (mark[iColumn] != stamp) {
            parallel = false;
            break;
          }
          const double small = newIsLarger ? element[j] : scatter[iColumn];
          const double large = newIsLarger ? scatter[iColumn] : element[j];
          if (std::fabs(large - up * small) > ROWRED_COEFTOL) {
            parallel = false;
            break;
          }
        }
        if (!parallel)
          continue;

        // The survivor is the row with the smaller coefficients, so the other
        // row's bounds are divided by |r| >= 1 on the way in. That keeps the
        // transferred bound from growing -- in particular a finite bound stays
        // finite and cannot cross the "infinite" threshold -- and, with |q|
        // == 1, leaves the lower row index as the survivor.
        const int survivor = newIsLarger ? jRow : iRow;
        const int victim = newIsLarger ? iRow : jRow;
        const double r = up; // row_victim == r * row_survivor, |r| >= 1

        double lo = workLower[victim];
        double hi = workUpper[victim];
        const bool loInfinite = lo <= -ROWRED_INFINITY;
        const bool hiInfinite = hi >= ROWRED_INFINITY;

        // Scaling by a negative r turns "at most hi" into "at least hi/r":
        // the two sides swap, and so does which of them is infinite. Getting
        // the *swap* wrong is the classic bug in this reduction, and it is
        // caught loudly by the tests. The explicit infinity *flags* are the
        // belt to that braces: because the survivor is always the smaller row,
        // |r| >= 1, so dividing a sentinel can only move it further from the
        // threshold and the arithmetic would in fact survive on its own here.
        // The flags mean the reduction stays correct if that survivor-choice
        // invariant is ever relaxed, which is the kind of change that would
        // otherwise break this silently.
        double transferLower, transferUpper;
        bool transferLowerInfinite, transferUpperInfinite;
        if (r > 0.0) {
          transferLower = lo / r;
          transferUpper = hi / r;
          transferLowerInfinite = loInfinite;
          transferUpperInfinite = hiInfinite;
        } else {
          transferLower = hi / r;
          transferUpper = lo / r;
          transferLowerInfinite = hiInfinite;
          transferUpperInfinite = loInfinite;
        }

        double newLower = workLower[survivor];
        double newUpper = workUpper[survivor];
        if (!transferLowerInfinite && transferLower > newLower)
          newLower = transferLower;
        if (!transferUpperInfinite && transferUpper < newUpper)
          newUpper = transferUpper;

        if (newLower > newUpper) {
          if (newLower - newUpper > ROWRED_FEASTOL
            && safeGreater(newLower, newUpper)) {
            if (logLevel >= 1)
              printf("  Row reductions: rows %d and %d are parallel (scale "
                     "%.12g) with incompatible bounds — infeasible.\n",
                survivor, victim, r);
            // Nothing is deleted or rewritten on the way out: every decision so
            // far is staged in `deleted` and `workLower`/`workUpper` and only
            // applied after both steps succeed, so returning here leaves the
            // model exactly as it arrived and the caller can still report on the
            // problem it handed in. Mutating this to delete the victim first is
            // caught by the tests; mutating it to *stage* a deletion is not, and
            // does not need to be -- a staged deletion that is never applied is
            // dead, which is the property this structure exists to give.
            infeasible_ = true;
            timeUsed_ = CoinWallclockTime() - t0;
            return false;
          }
          // Within tolerance of each other: snap to an exact equation rather
          // than leaving lower > upper by 1e-13, which trips asserts
          // downstream. Relaxing the lower bound (never the upper) keeps this
          // a widening, so no feasible point is lost.
          newLower = newUpper;
        }

        const bool tightened = std::fabs(newLower - workLower[survivor]) > ROWRED_EPSILON
          || std::fabs(newUpper - workUpper[survivor]) > ROWRED_EPSILON;
        if (tightened) {
          workLower[survivor] = newLower;
          workUpper[survivor] = newUpper;
          boundChanged[survivor] = 1;
          nBoundsTightened_++;
        }

        deleted[victim] = 1;
        const bool exactDuplicate = std::fabs(r - 1.0) <= ROWRED_EPSILON;
        if (exactDuplicate)
          nDuplicateRows_++;
        else
          nParallelRows_++;
        if (logLevel >= 3)
          printf("  Row reductions: row %d removed as %s of row %d "
                 "(scale %g)%s\n",
            victim, exactDuplicate ? "a duplicate" : "parallel to", survivor, r,
            tightened ? ", bounds tightened" : "");

        // If the new row displaced its representative, it takes its place so
        // that the rest of the bucket is still compared against something.
        if (victim == jRow)
          representative[rep] = iRow;
        merged = true;
      }

      if (!merged && numberRepresentatives < ROWRED_MAX_REPRESENTATIVES)
        representative[numberRepresentatives++] = iRow;
    }

    bucketStart = bucketEnd;
  }

  if (!nRowsRemoved() && !nBoundsTightened_) {
    timeUsed_ = CoinWallclockTime() - t0;
    return false;
  }

  // --- Self-check against a reference solution, when one is available -----
  //
  // Removing a row is a relaxation and cannot lose a feasible point; the risk
  // is entirely in the *tightened* bounds a merge leaves on a survivor. So
  // every row whose bounds changed is re-evaluated at the reference solution.
  //
  // Deliberately getRowCutDebuggerAlways() and NOT getRowCutDebugger(): the
  // latter gates on OsiRowCutDebugger::onOptimalPath(), which checks only
  // integer columns with a loose 1e-3 tolerance and has been observed to
  // return false while the current bounds are still fully consistent with the
  // reference solution -- see the matching comments in CbcModel.cpp,
  // CbcBoundPropagation.cpp and CbcCoefficientStrengthening.cpp. Nothing here
  // depends on which subtree we are in, so there is no reason to gate at all.
  // The debugSolution fallback covers "-debugCuts", which attaches its
  // debugger only after this phase has run.
  const OsiRowCutDebugger *debugger = solver->getRowCutDebuggerAlways();
  const double *optSol = debugger
    ? debugger->optimalSolution()
    : (debugSolution && debugNumberColumns == numberColumns ? debugSolution : nullptr);
  if (optSol) {
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (!boundChanged[iRow] || deleted[iRow])
        continue;
      double activity = 0.0;
      const CoinBigIndex start = rowStart[iRow];
      const int length = rowLength[iRow];
      for (CoinBigIndex j = start; j < start + length; j++)
        activity += element[j] * optSol[column[j]];
      const bool violatesUpper = workUpper[iRow] < ROWRED_INFINITY
        && activity > workUpper[iRow] + ROWRED_FEASTOL;
      const bool violatesLower = workLower[iRow] > -ROWRED_INFINITY
        && activity < workLower[iRow] - ROWRED_FEASTOL;
      if (violatesUpper || violatesLower) {
        printf("rowReductions BAD ROW: row %d activity=%.12g at reference "
               "solution but merged bounds are [%.12g,%.12g] (violates %s)\n",
          iRow, activity, workLower[iRow], workUpper[iRow],
          violatesUpper ? "upper" : "lower");
        fflush(stdout);
      }
    }
  }

  // --- Apply: bounds first (original indices), then one batched delete ----
  for (int iRow = 0; iRow < numberRows; iRow++) {
    if (!boundChanged[iRow] || deleted[iRow])
      continue;
    solver->setRowLower(iRow, workLower[iRow]);
    solver->setRowUpper(iRow, workUpper[iRow]);
  }

  const int totalRemoved = nRowsRemoved();
  if (totalRemoved) {
    std::vector< int > toRemove;
    toRemove.reserve(totalRemoved);
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (deleted[iRow])
        toRemove.push_back(iRow);
    }
    // Same mechanism CglCliqueStrengthening::removeDominatedRows() already
    // uses at this exact point of the pipeline: one batched deleteRows() on
    // pre-phase row indices. Safe here because nothing downstream has captured
    // them yet -- CbcSolver::preprocess() records truncateRows and
    // CglPreProcess builds originalRows() later, and the only earlier capture
    // is rowNames_, used for printing.
    solver->deleteRows(static_cast< int >(toRemove.size()), toRemove.data());
  }

  timeUsed_ = CoinWallclockTime() - t0;
  return true;
}
