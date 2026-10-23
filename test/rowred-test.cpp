/*
 * rowred-test — targeted tests for CbcRowReductions, the row-reduction step of
 * the pre-root-LP strengthening phase (step 4 of
 * CbcSolver::preRootLPStrenghtening(), per -rowReductions).
 *
 * The pass removes rows that cannot constrain the problem: rows all of whose
 * columns are fixed, and rows that are duplicates or scalar multiples of
 * another row -- with the survivor inheriting the intersection of the two rows'
 * bounds. So two things are worth asserting for every fixture:
 *
 *   1. the surviving rows, their count and their merged bounds are the ones
 *      worked out by hand; and
 *   2. no integer-feasible point is lost or gained -- checked by *enumerating*
 *      the whole box before and after, not by re-solving. A solver that reaches
 *      the same optimum tells us nothing about the points the reduced system
 *      silently stopped admitting, and for this pass "gained" matters as much
 *      as "lost": dropping a row is only sound if the survivor still says
 *      everything the dropped row said.
 *
 * Every fixture is small enough that Cbc would reach the right answer with the
 * pass disabled, so an objective-only test would pass on a build where the pass
 * did nothing. The bound assertions are therefore stated numerically, and the
 * two end-to-end fixtures are chosen so that the *specific* plausible bug
 * changes the optimum: a merge that keeps the looser of two duplicate bounds
 * reads 8 instead of 2, and a merge that forgets to swap the two bounds when
 * the scale factor is negative reads 6 instead of 8.
 *
 * The negative-scale cases carry most of the weight here. lambda = -1 is also
 * what exercises the hash's pivot tie-break (largest magnitude, ties broken by
 * smallest column index): a row and its negation have equal-magnitude candidate
 * pivots, and if they picked different ones they would normalise to different
 * vectors, never land in the same hash bucket, and never be compared -- so
 * every "1 row removed" assertion on a lambda = -1 fixture is also an assertion
 * that the tie-break works.
 *
 * Usage: rowred-test
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "CbcRowReductions.hpp"
#include "Cbc_C_Interface.h"
#include "CoinPackedMatrix.hpp"
#include "OsiClpSolverInterface.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool cond, const std::string &msg)
{
  if (!cond) {
    ++g_failures;
    printf("  FAIL: %s\n", msg.c_str());
  } else {
    printf("  ok:   %s\n", msg.c_str());
  }
}

void checkClose(double actual, double expected, double tol, const std::string &msg)
{
  double delta = fabs(actual - expected);
  if (delta > tol) {
    ++g_failures;
    printf("  FAIL: %s (actual=%.10g expected=%.10g delta=%.3e > tol=%.3e)\n",
      msg.c_str(), actual, expected, delta, tol);
  } else {
    printf("  ok:   %s (actual=%.10g expected=%.10g)\n",
      msg.c_str(), actual, expected);
  }
}

void checkEqInt(int actual, int expected, const std::string &msg)
{
  if (actual != expected) {
    ++g_failures;
    printf("  FAIL: %s (actual=%d expected=%d)\n", msg.c_str(), actual, expected);
  } else {
    printf("  ok:   %s (=%d)\n", msg.c_str(), actual);
  }
}

/* ---------------------------------------------------------------- fixtures */

/* A fixture stated the way it reads on paper: columns with bounds, integrality
 * and objective, then rows with explicit lower/upper bounds so that range rows,
 * one-sided rows and equations are equally expressible. (Same shape as
 * coefstr-test's Model, deliberately -- these two tests cover adjacent steps of
 * the same phase.) */
struct Model {
  std::vector< double > colLower, colUpper, obj;
  std::vector< bool > isInt;
  std::vector< std::vector< int > > rowIndex;
  std::vector< std::vector< double > > rowValue;
  std::vector< double > rowLower, rowUpper;

  int addCol(double lo, double up, bool integer, double c = 0.0)
  {
    colLower.push_back(lo);
    colUpper.push_back(up);
    isInt.push_back(integer);
    obj.push_back(c);
    return static_cast< int >(colLower.size()) - 1;
  }

  void addRow(const std::vector< int > &idx, const std::vector< double > &val,
    double lo, double up)
  {
    rowIndex.push_back(idx);
    rowValue.push_back(val);
    rowLower.push_back(lo);
    rowUpper.push_back(up);
  }

  void loadInto(OsiClpSolverInterface &si) const
  {
    const int nc = static_cast< int >(colLower.size());
    const int nr = static_cast< int >(rowIndex.size());
    CoinPackedMatrix m(false, 0, 0);
    m.setDimensions(0, nc);
    for (int i = 0; i < nr; i++)
      m.appendRow(static_cast< int >(rowIndex[i].size()),
        rowIndex[i].data(), rowValue[i].data());
    si.loadProblem(m, colLower.data(), colUpper.data(), obj.data(),
      rowLower.data(), rowUpper.data());
    for (int j = 0; j < nc; j++)
      if (isInt[j])
        si.setInteger(j);
    si.messageHandler()->setLogLevel(0);
  }
};

const double INF = 1.0e30;

/* Treat the solver's own sentinel the way the pass does. */
bool isInf(double v)
{
  return v >= 1.0e20;
}
bool isNegInf(double v)
{
  return v <= -1.0e20;
}

/* Assert a row's bounds, with "infinite" expressed as +/-INF. */
void checkRowBounds(OsiSolverInterface &si, int row, double lo, double up,
  const std::string &msg)
{
  if (row >= si.getNumRows()) {
    ++g_failures;
    printf("  FAIL: %s (row %d does not exist; only %d rows)\n",
      msg.c_str(), row, si.getNumRows());
    return;
  }
  const double actualLo = si.getRowLower()[row];
  const double actualUp = si.getRowUpper()[row];
  bool ok = true;
  if (isNegInf(lo))
    ok = ok && isNegInf(actualLo);
  else
    ok = ok && fabs(actualLo - lo) <= 1e-7;
  if (isInf(up))
    ok = ok && isInf(actualUp);
  else
    ok = ok && fabs(actualUp - up) <= 1e-7;
  if (!ok) {
    ++g_failures;
    printf("  FAIL: %s (row %d is [%.10g,%.10g], expected [%.10g,%.10g])\n",
      msg.c_str(), row, actualLo, actualUp, lo, up);
  } else {
    printf("  ok:   %s (row %d is [%.10g,%.10g])\n",
      msg.c_str(), row, actualLo, actualUp);
  }
}

/* Every integer point of the box that satisfies all rows, as a sorted list of
 * keys. Only usable on a tiny all-integer fixture -- which is the point: the
 * feasible sets before and after must be *identical*, and the only way to know
 * that is to look at all of them. */
std::vector< std::string > enumerateFeasible(OsiSolverInterface &si)
{
  const int nc = si.getNumCols();
  const int nr = si.getNumRows();
  const double *cl = si.getColLower();
  const double *cu = si.getColUpper();
  const double *rl = si.getRowLower();
  const double *ru = si.getRowUpper();
  const CoinPackedMatrix *m = si.getMatrixByRow();

  std::vector< std::string > out;
  std::vector< int > x(nc);
  for (int j = 0; j < nc; j++)
    x[j] = static_cast< int >(cl[j]);

  for (;;) {
    bool ok = true;
    for (int r = 0; r < nr && ok; r++) {
      double a = 0.0;
      const CoinBigIndex s = m->getVectorStarts()[r];
      const int len = m->getVectorLengths()[r];
      for (CoinBigIndex k = s; k < s + len; k++)
        a += m->getElements()[k] * x[m->getIndices()[k]];
      if (a > ru[r] + 1e-7 || a < rl[r] - 1e-7)
        ok = false;
    }
    if (ok) {
      std::string key;
      char buf[32];
      for (int j = 0; j < nc; j++) {
        snprintf(buf, sizeof(buf), "%d,", x[j]);
        key += buf;
      }
      out.push_back(key);
    }
    int j = 0; // odometer over the box
    for (; j < nc; j++) {
      if (x[j] < static_cast< int >(cu[j])) {
        x[j]++;
        break;
      }
      x[j] = static_cast< int >(cl[j]);
    }
    if (j == nc)
      break;
  }
  std::sort(out.begin(), out.end());
  return out;
}

/* Run the pass and assert the integer-feasible set is exactly preserved. */
void checkFeasibleSetPreserved(const Model &mo, const std::string &what)
{
  OsiClpSolverInterface before;
  mo.loadInto(before);
  const std::vector< std::string > setBefore = enumerateFeasible(before);

  OsiClpSolverInterface after;
  mo.loadInto(after);
  CbcRowReductions rr;
  rr.run(&after, after.messageHandler(), 0);
  const std::vector< std::string > setAfter = enumerateFeasible(after);

  if (setBefore == setAfter) {
    printf("  ok:   %s: integer-feasible set unchanged (%d points, %d rows -> "
           "%d rows)\n",
      what.c_str(), static_cast< int >(setBefore.size()), before.getNumRows(),
      after.getNumRows());
  } else {
    ++g_failures;
    printf("  FAIL: %s: integer-feasible set CHANGED (%d points -> %d points)\n",
      what.c_str(), static_cast< int >(setBefore.size()),
      static_cast< int >(setAfter.size()));
  }
}

/* ---------------------------------------------- 1. rows with all cols fixed */

/* x fixed at 2, z fixed at 1, y free to move.
 *
 *   row 0:  x + z <= 10        activity 3, can never bind  -> removed
 *   row 1:  x + z >= 1         activity 3, can never bind  -> removed
 *   row 2:  2 <= x + z <= 4    activity 3, inside the range -> removed
 *   row 3:  x + y <= 5         y is not fixed              -> kept
 *   row 4:  x + z = 3          activity 3, satisfied exactly -> removed
 */
void testAllFixedRows()
{
  printf("\n--- 1. rows with all columns fixed ---\n");

  Model mo;
  const int x = mo.addCol(2.0, 2.0, true);
  const int y = mo.addCol(0.0, 5.0, true);
  const int z = mo.addCol(1.0, 1.0, true);
  mo.addRow({ x, z }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, z }, { 1.0, 1.0 }, 1.0, INF);
  mo.addRow({ x, z }, { 1.0, 1.0 }, 2.0, 4.0);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
  mo.addRow({ x, z }, { 1.0, 1.0 }, 3.0, 3.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  const bool changed = rr.run(&si, si.messageHandler(), 0);

  check(changed, "the pass reports a change");
  check(!rr.isInfeasible(), "not reported infeasible");
  checkEqInt(rr.nFixedRows(), 4, "four all-fixed rows removed");
  checkEqInt(rr.nDuplicateRows() + rr.nParallelRows(), 0,
    "no duplicate/parallel removals (all-fixed rows are taken out first)");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -INF, 5.0, "the surviving row is x + y <= 5");

  checkFeasibleSetPreserved(mo, "all-fixed rows");
}

/* An all-fixed row whose constant activity violates its bounds proves the whole
 * model infeasible -- and nothing must be deleted once that is known, so the
 * caller can still report on the model it handed in. */
void testAllFixedRowInfeasible()
{
  printf("\n--- 2. all-fixed row that is violated -> infeasible ---\n");

  for (int side = 0; side < 2; side++) {
    Model mo;
    const int x = mo.addCol(2.0, 2.0, true);
    const int y = mo.addCol(0.0, 5.0, true);
    const int z = mo.addCol(1.0, 1.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
    if (side == 0) // activity 3 > 2
      mo.addRow({ x, z }, { 1.0, 1.0 }, -INF, 2.0);
    else // activity 3 < 8
      mo.addRow({ x, z }, { 1.0, 1.0 }, 8.0, INF);

    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    const bool changed = rr.run(&si, si.messageHandler(), 0);

    const std::string which = side == 0 ? "upper" : "lower";
    check(!changed, which + " bound violated: run() returns false");
    check(rr.isInfeasible(), which + " bound violated: isInfeasible() is true");
    checkEqInt(si.getNumRows(), 2,
      which + " bound violated: model left untouched (no rows deleted)");
  }
}

/* -------------------------------------------------- 3. exact duplicate rows */

/*   row 0:  x + y <= 10
 *   row 1:  x + y <= 6      identical support and coefficients (lambda = 1)
 *
 * Both rows pick x as their pivot (magnitudes tie at 1, tie broken by smallest
 * column index) and normalise to (1, 1), so they collide. Row 0 survives and
 * inherits min(10, 6) = 6.
 */
void testExactDuplicate()
{
  printf("\n--- 3. exact duplicate rows (lambda = 1) ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 6.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  const bool changed = rr.run(&si, si.messageHandler(), 0);

  check(changed, "the pass reports a change");
  checkEqInt(rr.nDuplicateRows(), 1, "one duplicate removed");
  checkEqInt(rr.nParallelRows(), 0, "counted as duplicate, not parallel");
  checkEqInt(rr.nBoundsTightened(), 1, "one row bound tightened");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -INF, 6.0, "survivor keeps the TIGHTER bound");

  checkFeasibleSetPreserved(mo, "exact duplicate");
}

/* The same pair with the columns listed in the opposite order inside the second
 * row. The per-nonzero hashes are combined order-independently precisely so
 * this still collides: CoinPackedMatrix makes no promise about the order of
 * indices within a vector, so an order-dependent combination would miss pairs
 * depending on how the model happened to be built. */
/* The same duplicate pair, but handed in with the column indices in different
 * orders. Note what this does and does not establish: OsiClp normalises each row
 * to ascending column order on the way in (verified directly -- a row built as
 * (2,0,1) reads back as (0,1,2) through both loadProblem and addRow), so by the
 * time the pass sees these two rows their storage order is already identical and
 * any hash combination at all would bucket them together. The fixture therefore
 * asserts that a caller may submit rows in any order, not that the combination
 * step is order-independent; that property is defence against a solver interface
 * that does not sort, and CoinPackedMatrix promises no order, so it cannot be
 * reached through OsiClp at all. Mutating the sum to an order-dependent
 * shift-mix accordingly leaves every fixture here passing.
 */
void testDuplicateUnsortedIndices()
{
  printf("\n--- 4. duplicate rows submitted with columns in different order ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  const int z = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y, z }, { 1.0, 2.0, 3.0 }, -INF, 12.0);
  mo.addRow({ z, x, y }, { 3.0, 1.0, 2.0 }, -INF, 9.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nDuplicateRows(), 1, "one duplicate removed despite the order");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -INF, 9.0, "survivor keeps the tighter bound");

  checkFeasibleSetPreserved(mo, "duplicate with reordered indices");
}

/* --------------------------------------------------- 5. parallel, lambda > 0 */

/*   row 0:  x + y <= 10
 *   row 1:  2x + 2y <= 14      = 2 * row 0
 *
 * Row 0 normalises to (1, 1) and so does row 1 (2/2, 2/2), so they collide.
 * The survivor is the row with the smaller coefficients -- row 0 -- and row 1's
 * bound comes across divided by 2: 14/2 = 7, tighter than 10.
 */
void testParallelPositive()
{
  printf("\n--- 5. parallel rows, lambda = 2 ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { 2.0, 2.0 }, -INF, 14.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 1, "one parallel row removed");
  checkEqInt(rr.nDuplicateRows(), 0, "counted as parallel, not duplicate");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -INF, 7.0, "survivor upper bound is 14/2 = 7");

  checkFeasibleSetPreserved(mo, "parallel lambda = 2");
}

/* Two range rows, so both bounds are finite and both have to be scaled, and
 * arranged so that the merge has to take one bound from each row. Every fixture
 * above happens to end up with the victim's bound winning, which an
 * implementation that simply *overwrites* the survivor's bounds with the
 * transferred ones would also satisfy; here that would be visibly wrong.
 *
 *   row 0:   1 <= x + y <= 9
 *   row 1:  -4 <= 2x + 2y <= 12    ->  -2 <= x + y <= 6
 *
 * Intersecting gives max(1, -2) = 1 from the SURVIVOR and min(9, 6) = 6 from the
 * victim.
 */
void testParallelRangeRowsIntersect()
{
  printf("\n--- 5b. parallel range rows: one bound from each ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 8.0, true);
  const int y = mo.addCol(0.0, 8.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, 1.0, 9.0);
  mo.addRow({ x, y }, { 2.0, 2.0 }, -4.0, 12.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 1, "one parallel row removed");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, 1.0, 6.0,
    "survivor keeps its OWN lower bound 1 and takes the victim's upper 12/2 = 6");

  checkFeasibleSetPreserved(mo, "parallel range rows");
}

/* --------------------------------------------------- 6. parallel, lambda < 0 */

/* The crux of the whole reduction, and the case most likely to be got wrong.
 *
 *   row 0:  x + y <= 10
 *   row 1:  -x - y <= -3       = -1 * row 0, i.e. x + y >= 3
 *
 * Scaling row 1 by lambda = -1 turns "at most -3" into "at least 3": the two
 * sides SWAP. So the survivor becomes 3 <= x + y <= 10.
 *
 * An implementation that forgets the swap computes upper = -3 / -1 = 3 and
 * leaves the lower bound at -infinity, producing x + y <= 3 -- which loses
 * every point with 4 <= x + y <= 10. The enumeration catches exactly that.
 *
 * This fixture is also the pivot tie-break test: row 1's two coefficients tie
 * at magnitude 1, as row 0's do, so the two rows only meet in the same hash
 * bucket because both resolve the tie to the smallest column index.
 */
void testParallelNegative()
{
  printf("\n--- 6. parallel rows, lambda = -1 (the bound swap) ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { -1.0, -1.0 }, -INF, -3.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 1,
    "one parallel row removed (so the pivot tie-break agreed on a column)");
  check(!rr.isInfeasible(), "not reported infeasible");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, 3.0, 10.0,
    "survivor is 3 <= x + y <= 10 -- the bounds SWAPPED on the way in");

  checkFeasibleSetPreserved(mo, "parallel lambda = -1");
}

/* lambda = -1 again, but arranged so that the pivot tie-break is the only thing
 * making the pair collide. The fixture above cannot show that: its rows are
 * (1, 1) and (-1, -1), which normalise to the same vector whichever of the two
 * tied entries is chosen as pivot, so it passes even with no tie-break at all.
 *
 *   row 0:  x - y <= 4         stored as (x: 1, y: -1)
 *   row 1: -x + y <= 2         stored as (y: 1, x: -1)   = -1 * row 0
 *
 * Both entries of both rows have magnitude 1, and the two rows list their
 * columns in opposite orders. Normalising row 1 by its x entry gives
 * (x: 1, y: -1) -- row 0 exactly. Normalising it by its y entry instead gives
 * (x: -1, y: 1), a different vector, a different hash and no collision. So only
 * a rule that resolves the tie by column index rather than by storage order
 * finds this pair.
 *
 * Merged: row 1 says x - y >= -2, so the survivor is -2 <= x - y <= 4.
 */
void testParallelNegativeTieBreak()
{
  printf("\n--- 6b. lambda = -1 where only the pivot tie-break collides ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 5.0, true);
  const int y = mo.addCol(0.0, 5.0, true);
  mo.addRow({ x, y }, { 1.0, -1.0 }, -INF, 4.0);
  mo.addRow({ y, x }, { 1.0, -1.0 }, -INF, 2.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 1,
    "one row removed -- the tie was broken by column index, not storage order");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -2.0, 4.0, "survivor is -2 <= x - y <= 4");

  checkFeasibleSetPreserved(mo, "lambda = -1 with reversed storage order");
}

/* Negative lambda where one side is genuinely infinite, so the *infinity* has
 * to swap sides too, not just the numbers.
 *
 *   row 0:  2x + 4y <= 20            [-inf, 20]
 *   row 1:  -x - 2y >= -30           [-30, +inf]   = -0.5 * row 0
 *
 * Row 1 has the smaller coefficients, so it survives and row 0's bounds come
 * across scaled by -2 with the sides swapped: row 0's finite upper bound 20
 * becomes a LOWER bound of 20 / -2 = -10, and row 0's infinite lower bound maps
 * to row 1's upper side, where it must stay infinite rather than being scaled
 * into a large finite number.
 *
 * Result: -10 <= -x - 2y <= +inf, which is exactly row 0 restated (x + 2y <= 10)
 * and strictly tighter than row 1's own -30.
 */
void testParallelNegativeWithInfinity()
{
  printf("\n--- 7. parallel rows, lambda < 0, one side infinite ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 2.0, 4.0 }, -INF, 20.0);
  mo.addRow({ x, y }, { -1.0, -2.0 }, -30.0, INF);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 1, "one parallel row removed");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -10.0, INF,
    "survivor is the SMALLER row, tightened to [-10, +inf]");

  checkFeasibleSetPreserved(mo, "parallel lambda = -0.5 with infinite side");
}

/* A large scale factor. The pass always divides the *larger* row's bounds by
 * |lambda| >= 1, so a finite bound can only shrink and an infinite one is
 * carried across as a flag rather than as arithmetic on the +/-1e20 sentinel.
 * That is what makes the sentinel convention safe here; this fixture locks the
 * large-|lambda| path in.
 *
 *   row 0:  1e6 x + 2e6 y <= 5e6      = 1e6 * row 1
 *   row 1:  x + 2y >= -1              [-1, +inf]
 *
 * Row 1 survives with upper bound 5e6 / 1e6 = 5, and its infinite side stays
 * infinite.
 */
void testParallelLargeRatio()
{
  printf("\n--- 8. parallel rows with a large scale factor ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0e6, 2.0e6 }, -INF, 5.0e6);
  mo.addRow({ x, y }, { 1.0, 2.0 }, -1.0, INF);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 1, "one parallel row removed");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, -1.0, 5.0,
    "survivor is the smaller row with upper bound 5e6/1e6 = 5");

  checkFeasibleSetPreserved(mo, "parallel with scale 1e6");
}

/* --------------------------------------------------- 9. contradictory merge */

/*   row 0:  x + y <= 3
 *   row 1:  -x - y <= -8       i.e. x + y >= 8
 *
 * Merging gives 8 <= x + y <= 3: infeasible, by a margin far outside both the
 * absolute and the relative tolerance. */
void testContradictoryMerge()
{
  printf("\n--- 9. parallel rows with incompatible bounds -> infeasible ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 3.0);
  mo.addRow({ x, y }, { -1.0, -1.0 }, -INF, -8.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  const bool changed = rr.run(&si, si.messageHandler(), 0);

  check(!changed, "run() returns false");
  check(rr.isInfeasible(), "isInfeasible() is true");
  checkEqInt(si.getNumRows(), 2, "model left untouched (no rows deleted)");
}

/* The same shape but with the two bounds crossing by less than the feasibility
 * tolerance. That is rounding noise, not infeasibility, so the merged bounds
 * snap to an exact equation instead of being left with lower > upper by 1e-7 --
 * which would trip asserts downstream. The snap relaxes the LOWER bound, so it
 * is a widening and cannot lose a feasible point.
 *
 *   row 0:  x + y <= 5
 *   row 1:  -x - y <= -5.0000001    i.e. x + y >= 5.0000001
 */
void testSnapToEquality()
{
  printf("\n--- 10. merged bounds cross within tolerance -> snap to equality ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
  mo.addRow({ x, y }, { -1.0, -1.0 }, -INF, -5.0000001);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  const bool changed = rr.run(&si, si.messageHandler(), 0);

  check(changed, "the pass reports a change");
  check(!rr.isInfeasible(), "NOT reported infeasible -- this is noise, not a proof");
  checkEqInt(si.getNumRows(), 1, "one row left");
  check(si.getRowLower()[0] <= si.getRowUpper()[0] + 1e-12,
    "survivor has lower <= upper (no hair-thin inverted range)");
  // Which way the snap goes matters, so this is asserted exactly rather than
  // within the usual 1e-7: relaxing the lower bound down to 5 is a widening and
  // keeps the six integer points with x + y = 5, whereas tightening the upper
  // bound up to 5.0000001 would give an equation no integer point satisfies.
  // A 1e-7 tolerance would accept both.
  checkClose(si.getRowUpper()[0], 5.0, 0.0,
    "the snap moved the LOWER bound: upper is still exactly 5");
  checkClose(si.getRowLower()[0], 5.0, 0.0,
    "the snap moved the LOWER bound: lower is now exactly 5, not 5.0000001");
  checkRowBounds(si, 0, 5.0, 5.0, "survivor snapped to the equation x + y = 5");
}

/* A column whose bounds are 1e-11 apart is *not* effectively fixed if its
 * coefficient is 1e9: the row's activity can still move by 1e-2, ten thousand
 * times the feasibility tolerance. So the all-fixed test has to weigh each near
 * fixing by its coefficient rather than looking at column bounds alone.
 *
 *   x in [0, 1e-11] continuous, coefficient 1e9 -- spread * |coef| = 1e-2
 *   row 0:  1e9 x <= 0.005
 *
 * The row genuinely restricts x, to <= 5e-12, which is half its upper bound.
 * But treating x as fixed at 0 gives a constant activity of 0, comfortably
 * inside the bound, and the row would be deleted -- losing that restriction
 * entirely. The weighted test keeps it.
 *
 * The second half guards the other side of the same test: a movement that really
 * is negligible must still count as fixed, so the test cannot be passed by
 * simply never treating anything as fixed.
 *
 * (The spreads are stated as [0, eps] rather than [1, 1+eps] on purpose --
 * 1.0 + 1e-11 - 1.0 rounds to 1.0000000827e-11, so a bound pair straddling 1
 * cannot express a spread of exactly 1e-11 and the fixture would be testing the
 * rounding rather than the rule.)
 */
void testNearlyFixedColumnWithLargeCoefficient()
{
  printf("\n--- 10b. near-fixed column is weighed by its coefficient ---\n");

  {
    Model mo;
    const int x = mo.addCol(0.0, 1.0e-11, false);
    mo.addRow({ x }, { 1.0e9 }, -INF, 0.005);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    const bool changed = rr.run(&si, si.messageHandler(), 0);
    check(!changed, "spread 1e-11 with coefficient 1e9: no change reported");
    checkEqInt(rr.nFixedRows(), 0, "not counted as an all-fixed row");
    checkEqInt(si.getNumRows(), 1, "the row is kept");
  }

  { // the same shape with a coefficient of 1: genuinely negligible, so removed
    Model mo;
    const int x = mo.addCol(0.0, 1.0e-11, false);
    mo.addRow({ x }, { 1.0 }, -INF, 0.5);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(rr.nFixedRows(), 1, "spread 1e-11 with coefficient 1: IS fixed");
    checkEqInt(si.getNumRows(), 0, "the row is removed");
  }
}

/* An all-fixed row whose activity misses its bound by less than the feasibility
 * tolerance sits in the sliver between "clearly violated" (infeasible) and
 * "comfortably satisfied" (removable), and must be left in place. Deleting it
 * would be unsound in the direction that matters: the row is violated, so
 * dropping it turns an infeasible model into a feasible-looking one.
 *
 *   x fixed at 2, z fixed at 1, activity = 3
 *   row 0:  x + z <= 3 - 1e-7      violated by 1e-7 < feastol  -> kept
 *   row 1:  x + z >= 3 + 1e-7      violated by 1e-7 < feastol  -> kept
 *
 * Neither is reported infeasible either -- 1e-7 is inside the tolerance, so the
 * pass declines to make that call and hands the decision to the LP, which is the
 * asymmetry the pass is built around (a wrong redundancy claim drops one
 * constraint; a wrong infeasibility claim throws away the whole model).
 *
 * This is the fixture that makes the satisfaction test load-bearing: without it,
 * replacing `if (safeUpper && safeLower)` with `if (true)` passes everything,
 * because every other all-fixed fixture is either comfortably satisfied or
 * violated by enough to be caught by the infeasibility branch first.
 */
void testAllFixedRowMarginallyViolated()
{
  printf("\n--- 10c. all-fixed row violated by less than feastol is kept ---\n");

  for (int side = 0; side < 2; side++) {
    Model mo;
    const int x = mo.addCol(2.0, 2.0, true);
    const int z = mo.addCol(1.0, 1.0, true);
    if (side == 0)
      mo.addRow({ x, z }, { 1.0, 1.0 }, -INF, 3.0 - 1.0e-7);
    else
      mo.addRow({ x, z }, { 1.0, 1.0 }, 3.0 + 1.0e-7, INF);

    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    const bool changed = rr.run(&si, si.messageHandler(), 0);

    const std::string which = side == 0 ? "upper" : "lower";
    check(!changed, which + " missed by 1e-7: no change reported");
    check(!rr.isInfeasible(),
      which + " missed by 1e-7: NOT declared infeasible (inside tolerance)");
    checkEqInt(rr.nFixedRows(), 0, which + " missed by 1e-7: not counted");
    checkEqInt(si.getNumRows(), 1, which + " missed by 1e-7: the row is kept");
  }
}

/* ------------------------------------------- 11. three mutually parallel rows */

/* All three normalise to (1, 1) and land in one bucket, so the whole class
 * collapses to a single row carrying every bound:
 *
 *   row 0:  x + y <= 10
 *   row 1:  2x + 2y <= 14      ->  x + y <= 7
 *   row 2:  -x - y <= -3       ->  x + y >= 3
 *
 * leaving 3 <= x + y <= 7.
 */
void testThreeMutuallyParallel()
{
  printf("\n--- 11. a whole parallel class collapses to one row ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { 2.0, 2.0 }, -INF, 14.0);
  mo.addRow({ x, y }, { -1.0, -1.0 }, -INF, -3.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nParallelRows(), 2, "two of the three rows removed");
  checkEqInt(si.getNumRows(), 1, "one row left");
  checkRowBounds(si, 0, 3.0, 7.0, "survivor carries every bound: 3 <= x + y <= 7");

  checkFeasibleSetPreserved(mo, "three mutually parallel rows");
}

/* ------------------------------------------------------- 12. must not fire */

/* Eight rows the reduction must decline, each for its own reason. Any one of
 * them firing would remove a constraint the model still needs. */
void testMustNotFire()
{
  printf("\n--- 12. rows the reduction must decline ---\n");

  { // (a) singleton rows: a bound in disguise, and bound propagation's job
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x }, { 1.0 }, -INF, 3.0);
    mo.addRow({ x }, { 1.0 }, -INF, 5.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(si.getNumRows(), 2, "(a) two singleton rows both kept");
  }

  { // (b) same support, coefficients not proportional
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
    mo.addRow({ x, y }, { 1.0, 2.0 }, -INF, 5.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(si.getNumRows(), 2, "(b) same support but not proportional: both kept");
  }

  { // (c) proportional to within 1e-7 only -- far outside the 1e-9 coefficient
    // tolerance. The hash may well bucket these together; the explicit
    // coefficient-by-coefficient verification is what rejects them.
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
    mo.addRow({ x, y }, { 1.0, 1.0000001 }, -INF, 4.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(si.getNumRows(), 2,
      "(c) proportional only to 1e-7: both kept (verification rejects)");
  }

  { // (d) different support of the same size
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    const int z = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
    mo.addRow({ x, z }, { 1.0, 1.0 }, -INF, 4.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(si.getNumRows(), 2, "(d) different support: both kept");
  }

  { // (e) one row's support is a strict subset of the other's
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    const int z = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
    mo.addRow({ x, y, z }, { 1.0, 1.0, 1.0 }, -INF, 4.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(si.getNumRows(), 2, "(e) subset support: both kept");
  }

  { // (f) partly-fixed row: deliberately out of scope, folding the fixed
    // column into the bounds is bound propagation's and CglPreProcess' job
    Model mo;
    const int x = mo.addCol(2.0, 2.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 5.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    const bool changed = rr.run(&si, si.messageHandler(), 0);
    check(!changed, "(f) partly-fixed row: no change reported");
    checkEqInt(si.getNumRows(), 1, "(f) partly-fixed row kept");
  }

  { // (g) a model with nothing to find at all: the pass must report no change,
    // which is the path it spends nearly all its time on in practice
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    const int z = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 2.0 }, -INF, 5.0);
    mo.addRow({ y, z }, { 3.0, 1.0 }, -INF, 4.0);
    mo.addRow({ x, z }, { 1.0, -7.0 }, 1.0, 9.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    const bool changed = rr.run(&si, si.messageHandler(), 0);
    check(!changed, "(g) nothing to find: no change reported");
    checkEqInt(rr.nRowsRemoved(), 0, "(g) nothing to find: no rows removed");
    checkEqInt(si.getNumRows(), 3, "(g) all three rows kept");
  }

  { // (h) an empty model must not crash or claim anything
    OsiClpSolverInterface si;
    Model mo;
    mo.addCol(0.0, 1.0, true);
    mo.loadInto(si);
    CbcRowReductions rr;
    const bool changed = rr.run(&si, si.messageHandler(), 0);
    check(!changed && !rr.isInfeasible(), "(h) model with no rows: no change, feasible");
  }
}

/* ------------------------------------------------------ 13. duplicate equations */

/* Two identical equations, and two equations related by lambda = -1. Both
 * collapse, and the survivor stays an equation. The lambda = -1 case is the one
 * where a missing swap is invisible in the bounds (an equation scaled by -1 is
 * still an equation) but the *numbers* still have to move: -x - y = -4 must
 * become x + y = 4, not x + y = -4. */
void testDuplicateEquations()
{
  printf("\n--- 13. duplicate and negated equations ---\n");

  { // identical
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, 4.0, 4.0);
    mo.addRow({ x, y }, { 1.0, 1.0 }, 4.0, 4.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(rr.nDuplicateRows(), 1, "identical equations: one removed");
    checkEqInt(si.getNumRows(), 1, "identical equations: one row left");
    checkRowBounds(si, 0, 4.0, 4.0, "identical equations: survivor is x + y = 4");
    checkFeasibleSetPreserved(mo, "identical equations");
  }

  { // negated: x + y = 4 and -x - y = -4 are the same constraint
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, 4.0, 4.0);
    mo.addRow({ x, y }, { -1.0, -1.0 }, -4.0, -4.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(rr.nParallelRows(), 1, "negated equations: one removed");
    check(!rr.isInfeasible(), "negated equations: NOT infeasible");
    checkEqInt(si.getNumRows(), 1, "negated equations: one row left");
    checkRowBounds(si, 0, 4.0, 4.0, "negated equations: survivor is still x + y = 4");
    checkFeasibleSetPreserved(mo, "negated equations");
  }

  { // an equation and a matching inequality: the equation wins outright
    Model mo;
    const int x = mo.addCol(0.0, 10.0, true);
    const int y = mo.addCol(0.0, 10.0, true);
    mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 9.0);
    mo.addRow({ x, y }, { 2.0, 2.0 }, 8.0, 8.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcRowReductions rr;
    rr.run(&si, si.messageHandler(), 0);
    checkEqInt(si.getNumRows(), 1, "equation vs inequality: one row left");
    checkRowBounds(si, 0, 4.0, 4.0,
      "equation vs inequality: survivor is the equation x + y = 4");
    checkFeasibleSetPreserved(mo, "equation absorbing a weaker inequality");
  }
}

/* ------------------------------------------ 14. order independence of the result */

/* The same three constraints presented in reverse row order. The surviving
 * system must describe the same feasible region -- which is the property that
 * matters, and a stronger claim than "the same row index survived". */
void testRowOrderIndependence()
{
  printf("\n--- 14. result does not depend on row order ---\n");

  Model forward;
  const int x1 = forward.addCol(0.0, 10.0, true);
  const int y1 = forward.addCol(0.0, 10.0, true);
  forward.addRow({ x1, y1 }, { 1.0, 1.0 }, -INF, 10.0);
  forward.addRow({ x1, y1 }, { 2.0, 2.0 }, -INF, 14.0);
  forward.addRow({ x1, y1 }, { -1.0, -1.0 }, -INF, -3.0);

  Model reverse;
  const int x2 = reverse.addCol(0.0, 10.0, true);
  const int y2 = reverse.addCol(0.0, 10.0, true);
  reverse.addRow({ x2, y2 }, { -1.0, -1.0 }, -INF, -3.0);
  reverse.addRow({ x2, y2 }, { 2.0, 2.0 }, -INF, 14.0);
  reverse.addRow({ x2, y2 }, { 1.0, 1.0 }, -INF, 10.0);

  OsiClpSolverInterface a, b;
  forward.loadInto(a);
  reverse.loadInto(b);
  CbcRowReductions rrA, rrB;
  rrA.run(&a, a.messageHandler(), 0);
  rrB.run(&b, b.messageHandler(), 0);

  checkEqInt(a.getNumRows(), b.getNumRows(),
    "same number of rows survive either way");
  checkEqInt(rrA.nRowsRemoved(), rrB.nRowsRemoved(),
    "same number of rows removed either way");
  check(enumerateFeasible(a) == enumerateFeasible(b),
    "the two reduced systems describe the same feasible region");

  checkFeasibleSetPreserved(reverse, "reversed row order");
}

/* --------------------------------------------------------- 15. idempotence */

/* A second pass over an already-reduced model must find nothing. If it does,
 * the first pass left a reduction on the table for no reason, or -- worse --
 * the merge produced a new pair it should have handled in one go. */
void testIdempotent()
{
  printf("\n--- 15. a second pass finds nothing ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  const int z = mo.addCol(2.0, 2.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { 2.0, 2.0 }, -INF, 14.0);
  mo.addRow({ x, y }, { -1.0, -1.0 }, -INF, -3.0);
  mo.addRow({ z }, { 1.0 }, -INF, 5.0);
  mo.addRow({ x, z }, { 1.0, 1.0 }, -INF, 9.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions first;
  check(first.run(&si, si.messageHandler(), 0), "first pass finds something");
  const int rowsAfterFirst = si.getNumRows();

  CbcRowReductions second;
  const bool changedAgain = second.run(&si, si.messageHandler(), 0);
  check(!changedAgain, "second pass finds nothing");
  checkEqInt(second.nRowsRemoved(), 0, "second pass removes no rows");
  checkEqInt(si.getNumRows(), rowsAfterFirst, "row count unchanged by the second pass");
}

/* ------------------------------------------------- 16. accounting and timing */

void testAccounting()
{
  printf("\n--- 16. counters add up ---\n");

  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  const int z = mo.addCol(3.0, 3.0, true);
  mo.addRow({ z }, { 2.0 }, -INF, 100.0); // all-fixed (singleton, but Step A
                                          // runs before the length filter)
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 8.0); // duplicate
  mo.addRow({ x, y }, { 3.0, 3.0 }, -INF, 21.0); // parallel, lambda = 3
  mo.addRow({ x, z }, { 1.0, 1.0 }, -INF, 9.0); // partly fixed: kept

  OsiClpSolverInterface si;
  mo.loadInto(si);
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);

  checkEqInt(rr.nFixedRows(), 1, "one all-fixed row");
  checkEqInt(rr.nDuplicateRows(), 1, "one duplicate row");
  checkEqInt(rr.nParallelRows(), 1, "one parallel row");
  checkEqInt(rr.nRowsRemoved(),
    rr.nFixedRows() + rr.nDuplicateRows() + rr.nParallelRows(),
    "nRowsRemoved() is the sum of the three");
  checkEqInt(si.getNumRows(), 5 - rr.nRowsRemoved(),
    "row count dropped by exactly nRowsRemoved()");
  check(rr.timeUsed() >= 0.0 && rr.timeUsed() < 5.0,
    "timeUsed() is a sane wall-clock figure");
  // min(10, 8, 21/3 = 7) = 7
  checkRowBounds(si, 0, -INF, 7.0, "survivor carries the tightest of the three");

  checkFeasibleSetPreserved(mo, "mixed fixture");
}

/* ------------------------------------------------------- 17. end to end */

/* Through the CLI parameter, on two fixtures chosen so that the plausible bug
 * changes the reported optimum -- otherwise the assertion would pass on a build
 * where the pass did nothing at all.
 *
 * Fixture A (duplicate rows, tighter bound on the second):
 *     max y  s.t.  x + y <= 10,  x + y <= 4,  x = 2,  y in [0,10] integer
 *   Correct: y = 2. A merge that keeps the LOOSER bound reads y = 8.
 *
 * Fixture B (lambda = -1):
 *     max y  s.t.  x + y <= 10,  -x - y <= -8,  x = 2,  y in [0,10] integer
 *   Correct: 8 <= x + y <= 10 with x = 2, so y = 8. A merge that forgets the
 *   negative-lambda swap turns the pair into x + y <= 8 and reads y = 6.
 *
 * Both are also run with -rowReductions off, which must give the same answer:
 * the reduction changes the formulation, never the optimum.
 */
void testParameterEndToEnd()
{
  printf("\n--- 17. end to end: -rowReductions on/off ---\n");

  for (int fixture = 0; fixture < 2; fixture++) {
    const char *name = fixture == 0 ? "duplicate rows" : "lambda = -1 rows";
    const double expected = fixture == 0 ? 2.0 : 8.0;
    for (int off = 0; off < 2; off++) {
      Cbc_Model *model = Cbc_newModel();
      Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 0);
      if (off)
        Cbc_setParameter(model, "rowReductions", "off");
      Cbc_setParameter(model, "preprocess", "off");

      Cbc_addCol(model, "x", 0.0, 10.0, 0.0, 1, 0, NULL, NULL);
      Cbc_addCol(model, "y", 0.0, 10.0, -1.0, 1, 0, NULL, NULL); // min -y = max y
      {
        int cols[2] = { 0, 1 };
        double coefs[2] = { 1.0, 1.0 };
        Cbc_addRow(model, "r0", 2, cols, coefs, 'L', 10.0);
      }
      if (fixture == 0) {
        int cols[2] = { 0, 1 };
        double coefs[2] = { 1.0, 1.0 };
        Cbc_addRow(model, "r1", 2, cols, coefs, 'L', 4.0);
      } else {
        int cols[2] = { 0, 1 };
        double coefs[2] = { -1.0, -1.0 };
        Cbc_addRow(model, "r1", 2, cols, coefs, 'L', -8.0);
      }
      { // x = 2
        int cols[1] = { 0 };
        double coefs[1] = { 1.0 };
        Cbc_addRow(model, "pin", 1, cols, coefs, 'E', 2.0);
      }

      Cbc_solve(model);
      const std::string what = std::string(name) + ": -rowReductions "
        + (off ? "off" : "on") + " gives max y = "
        + (fixture == 0 ? "2" : "8");
      checkClose(-Cbc_getObjValue(model), expected, 1e-6, what);
      Cbc_deleteModel(model);
    }
  }

  // And the point of the pass, asserted directly: the model handed to the root
  // LP really is smaller.
  Model mo;
  const int x = mo.addCol(0.0, 10.0, true);
  const int y = mo.addCol(0.0, 10.0, true);
  mo.addRow({ x, y }, { 1.0, 1.0 }, -INF, 10.0);
  mo.addRow({ x, y }, { 2.0, 2.0 }, -INF, 14.0);
  mo.addRow({ x, y }, { -1.0, -1.0 }, -INF, -3.0);
  OsiClpSolverInterface si;
  mo.loadInto(si);
  const int rowsBefore = si.getNumRows();
  CbcRowReductions rr;
  rr.run(&si, si.messageHandler(), 0);
  check(si.getNumRows() < rowsBefore,
    "the model the pass hands on has strictly fewer rows");
  printf("        (rows: %d before, %d after)\n", rowsBefore, si.getNumRows());
}

} // end anonymous namespace

int main(int argc, char **argv)
{
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      printf("rowred-test — tests for CbcRowReductions (removal of all-fixed,\n"
             "duplicate and parallel rows, step 4 of the pre-root-LP\n"
             "strengthening phase).\n\n"
             "Usage: rowred-test\n\n"
             "Takes no options and needs no data files: every fixture and its\n"
             "expected surviving rows and merged bounds are built in-process\n"
             "and worked out by hand in the comment above each test.\n\n"
             "Exit code: 0 = all checks passed, 1 = one or more failed.\n");
      return 0;
    }
  }

  printf("=== rowred-test: row reductions (all-fixed, duplicate, parallel) ===\n");

  testAllFixedRows();
  testAllFixedRowInfeasible();
  testExactDuplicate();
  testDuplicateUnsortedIndices();
  testParallelPositive();
  testParallelRangeRowsIntersect();
  testParallelNegative();
  testParallelNegativeTieBreak();
  testParallelNegativeWithInfinity();
  testParallelLargeRatio();
  testContradictoryMerge();
  testSnapToEquality();
  testNearlyFixedColumnWithLargeCoefficient();
  testAllFixedRowMarginallyViolated();
  testThreeMutuallyParallel();
  testMustNotFire();
  testDuplicateEquations();
  testRowOrderIndependence();
  testIdempotent();
  testAccounting();
  testParameterEndToEnd();

  printf("\n");
  if (g_failures) {
    printf("RESULT: %d check(s) FAILED\n", g_failures);
    return 1;
  }
  printf("RESULT: all checks passed\n");
  return 0;
}
