/*
 * coefstr-test — targeted tests for CbcCoefficientStrengthening, the
 * coefficient tightening ("big-M strengthening") step of the pre-root-LP
 * strengthening phase (step 3 of CbcSolver::preRootLPStrenghtening()).
 *
 * The pass shrinks an integer coefficient that exceeds its row's slack down to
 * that slack, adjusting the right-hand side to compensate. That tightens the LP
 * relaxation while leaving the integer-feasible set alone, so two things are
 * worth asserting for every fixture:
 *
 *   1. the new coefficient and right-hand side are the ones worked out by hand,
 *      and the LP relaxation bound moves by the predicted amount; and
 *   2. no integer-feasible point is lost -- checked by *enumerating* the whole
 *      box, not by re-solving. A solver that reaches the same optimum tells us
 *      nothing about the points it silently stopped admitting.
 *
 * Every fixture here is small enough that Cbc would reach the right answer with
 * the pass disabled, so an objective-only test would pass on a build where the
 * pass did nothing. Each assertion was therefore checked by mutation: the
 * coefficient and rhs values, the "must not fire" guards, and the enumeration
 * are each strong enough to fail on a plausible wrong implementation -- rounding
 * the slack down, adjusting the rhs off the wrong bound, firing on range rows,
 * treating a continuous column as integer, or leaving the solver's cached row
 * copy stale.
 *
 * Usage: coefstr-test
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "CbcCoefficientStrengthening.hpp"
#include "Cbc_C_Interface.h"
#include "CoinPackedMatrix.hpp"
#include "OsiClpSolverInterface.hpp"

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
 * and objective, then rows with explicit lower/upper bounds so that range rows
 * and one-sided rows are equally expressible. */
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

/* Coefficient of (row, col) read through the solver's ROW view. Reading it this
 * way is deliberate: the cached row copy is exactly what a
 * modifyCoefficient()-based implementation would leave stale. */
double rowCoef(OsiSolverInterface &si, int row, int col)
{
  const CoinPackedMatrix *m = si.getMatrixByRow();
  const CoinBigIndex start = m->getVectorStarts()[row];
  const int length = m->getVectorLengths()[row];
  const int *idx = m->getIndices();
  const double *el = m->getElements();
  for (CoinBigIndex k = start; k < start + length; k++)
    if (idx[k] == col)
      return el[k];
  return 0.0;
}

/* LP relaxation objective: continuous relaxation only, no cuts, no branching. */
double lpBound(const OsiClpSolverInterface &si)
{
  OsiClpSolverInterface copy(si);
  copy.messageHandler()->setLogLevel(0);
  copy.initialSolve();
  return copy.isProvenOptimal() ? copy.getObjValue() : COIN_DBL_MAX;
}

/* Every integer point of the box that satisfies all rows, as a list of keys.
 * Only usable on a tiny all-integer fixture -- which is the point: the feasible
 * sets before and after must be *identical*, and the only way to know that is
 * to look at all of them. */
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
  return out;
}

/* ------------------------------------------------------- 1. textbook big-M */

/* The canonical fixed-charge row. y is the on/off switch for x:
 *
 *     x - 100 y <= 0,   0 <= x <= 7 continuous,  y in {0,1}
 *
 * Max activity is 7 (x at 7, y at its lower bound 0 since its coefficient is
 * negative), so slack = 7 - 0 = 7. |-100| > 7, so the coefficient shrinks to
 * -7, and the rhs is adjusted off y's LOWER bound: rhs -= (-100 + 7) * 0 = 0.
 *
 * The LP relaxation is the whole point. The row forces y >= x/100, so
 * *minimising* y with x pinned at 7 by a second row gives 0.07 originally and 1
 * after tightening -- the value integrality implies anyway, but now visible to
 * the very first LP rather than only after cuts.
 */
void testTextbookBigM()
{
  printf("\n--- 1. textbook big-M row: x - 100 y <= 0 ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 7.0, false, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0); // minimise y
  mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
  mo.addRow({ x }, { 1.0 }, 7.0, 7.0); // pin x at 7 so the big-M row bites

  OsiClpSolverInterface si;
  mo.loadInto(si);
  const double before = lpBound(si);

  CbcCoefficientStrengthening cs;
  const bool changed = cs.run(&si, si.messageHandler(), 0);

  check(changed, "run() reports a change");
  checkEqInt(cs.nCoefficients(), 1, "one coefficient tightened");
  checkEqInt(cs.nRows(), 1, "one row touched");
  checkClose(rowCoef(si, 0, y), -7.0, 1e-9, "coefficient of y: -100 -> -7 (slack)");
  checkClose(rowCoef(si, 0, x), 1.0, 1e-9, "coefficient of x unchanged (1 <= slack)");
  checkClose(si.getRowUpper()[0], 0.0, 1e-9, "rhs unchanged (y's lower bound is 0)");
  checkClose(before, 0.07, 1e-6, "LP bound before: y only forced up to 7/100");
  checkClose(lpBound(si), 1.0, 1e-6, "LP bound after: y forced up to 1");
  check(lpBound(si) > before + 0.5, "LP bound is strictly, substantially tighter");
}

/* -------------------------------------- 2. general integer, not just binary */

/* The case Cgl's own tightener skips: its gate is `difference <= 1.0`, i.e.
 * binary columns only, so a general integer keeps its oversized coefficient.
 * Generalising past that gate is precisely what this pass adds.
 *
 *     30 z + 4 v <= 42,   z integer in [0,2],  v integer in [0,2]
 *
 * Max activity = 60 + 8 = 68, slack = 68 - 42 = 26.
 *   - z: |30| > 26, so it shrinks to 26 and the rhs moves off z's UPPER bound
 *     (positive coefficient): rhs = 42 - (30 - 26) * 2 = 34.
 *   - v: |4| <= 26, untouched.
 *
 * The LP cap on z tightens from 42/30 = 1.4 to 34/26 ~= 1.3077. Note the rule
 * needs |a| > slack *strictly*: widening v's bound to 3 gives slack 30 and
 * changes nothing.
 */
void testGeneralInteger()
{
  printf("\n--- 2. general integer (not binary): 30 z + 4 v <= 42 ---\n");
  Model mo;
  const int z = mo.addCol(0.0, 2.0, true, -1.0); // maximise z
  const int v = mo.addCol(0.0, 2.0, true, 0.0);
  mo.addRow({ z, v }, { 30.0, 4.0 }, -COIN_DBL_MAX, 42.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  const std::vector< std::string > feasibleBefore = enumerateFeasible(si);
  const double before = lpBound(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);

  checkEqInt(cs.nCoefficients(), 1, "one coefficient tightened (general integer z)");
  checkClose(rowCoef(si, 0, z), 26.0, 1e-9, "coefficient of z: 30 -> 26 (slack)");
  checkClose(rowCoef(si, 0, v), 4.0, 1e-9, "coefficient of v unchanged (4 <= slack)");
  checkClose(si.getRowUpper()[0], 34.0, 1e-9,
    "rhs adjusted off z's UPPER bound: 42 - (30-26)*2 = 34");
  checkClose(before, -1.4, 1e-6, "LP bound before: z <= 42/30 = 1.4");
  checkClose(lpBound(si), -34.0 / 26.0, 1e-6, "LP bound after: z <= 34/26");
  check(lpBound(si) > before, "LP bound is strictly tighter");

  check(enumerateFeasible(si) == feasibleBefore,
    "integer-feasible set unchanged (full enumeration)");
  check(!feasibleBefore.empty(), "enumeration found feasible points at all");
}

/* ------------------------------------------------- 3. negative coefficient */

/* A negative oversized coefficient must adjust the rhs off the column's LOWER
 * bound, since that is where maximum activity puts it. Getting this backwards
 * is the easy mistake, and it cuts off feasible points.
 *
 *     4 v - 30 z <= -26,   v integer in [0,2],  z integer in [1,3]
 *
 * Max activity: v at 2 (+8), z at its lower bound 1 (-30) => -22.
 * slack = -22 - (-26) = 4, and |-30| > 4, so -30 shrinks to -4 with
 * rhs -= (-30 - (-4)) * 1 = +26  =>  rhs = -26 + 26 = 0.
 */
void testNegativeCoefficient()
{
  printf("\n--- 3. negative coefficient: 4 v - 30 z <= -26 ---\n");
  Model mo;
  const int v = mo.addCol(0.0, 2.0, true, 0.0);
  const int z = mo.addCol(1.0, 3.0, true, 1.0); // minimise z
  mo.addRow({ v, z }, { 4.0, -30.0 }, -COIN_DBL_MAX, -26.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  const std::vector< std::string > feasibleBefore = enumerateFeasible(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);

  checkEqInt(cs.nCoefficients(), 1, "one coefficient tightened");
  checkClose(rowCoef(si, 0, z), -4.0, 1e-9, "coefficient of z: -30 -> -4");
  checkClose(si.getRowUpper()[0], 0.0, 1e-9,
    "rhs adjusted off z's LOWER bound: -26 - (-30+4)*1 = 0");
  check(enumerateFeasible(si) == feasibleBefore,
    "integer-feasible set unchanged (full enumeration)");
  check(!feasibleBefore.empty(), "enumeration found feasible points at all");
}

/* ------------------------------------------------------------- 4. >= rows */

/* A >= row is the same rule after negating the whole row (scale = -1), so the
 * mirror of test 1 must come back with the mirrored coefficient:
 *
 *     -x + 100 y >= 0   negates to   x - 100 y <= 0
 *
 * y's coefficient must therefore become +7, and the row's LOWER bound (not its
 * upper) is the one that may move.
 */
void testGreaterEqualRow()
{
  printf("\n--- 4. >= row (the scale = -1 path): -x + 100 y >= 0 ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 7.0, false, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0); // minimise y
  mo.addRow({ x, y }, { -1.0, 100.0 }, 0.0, COIN_DBL_MAX);
  mo.addRow({ x }, { 1.0 }, 7.0, 7.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  const double before = lpBound(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);

  checkEqInt(cs.nCoefficients(), 1, "one coefficient tightened on the >= row");
  checkClose(rowCoef(si, 0, y), 7.0, 1e-9, "coefficient of y: 100 -> 7");
  checkClose(si.getRowLower()[0], 0.0, 1e-9, "row lower bound unchanged");
  check(si.getRowUpper()[0] > 1.0e20, "row is still one-sided (>=)");
  checkClose(before, 0.07, 1e-6, "LP bound before: y only forced up to 0.07");
  checkClose(lpBound(si), 1.0, 1e-6, "LP bound after: y forced up to 1");
}

/* --------------------------------------------------------- 5. must NOT fire */

/* Each fixture below has an oversized coefficient that the rule must
 * nevertheless leave alone, and each guard is a real hazard: a range row has
 * two slacks and neither justifies the change alone; an infinite max activity
 * means there is no slack at all; a singleton row's "tightening" is a bound
 * change, which is bound propagation's job; a continuous column cannot be
 * stepped by 1, which is what the redundancy argument needs; and a fixed or
 * sub-unit-range column has no step either.
 */
void testMustNotFire()
{
  printf("\n--- 5. must not fire ---\n");

  { // range row: -1 <= x - 100 y <= 0
    Model mo;
    const int x = mo.addCol(0.0, 7.0, false);
    const int y = mo.addCol(0.0, 1.0, true);
    mo.addRow({ x, y }, { 1.0, -100.0 }, -1.0, 0.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "range row: no change");
    checkClose(rowCoef(si, 0, y), -100.0, 1e-9, "range row: coefficient intact");
  }

  { // equality row: also two finite sides
    Model mo;
    const int x = mo.addCol(0.0, 7.0, false);
    const int y = mo.addCol(0.0, 1.0, true);
    mo.addRow({ x, y }, { 1.0, -100.0 }, 0.0, 0.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "equality row: no change");
    checkClose(rowCoef(si, 0, y), -100.0, 1e-9, "equality row: coefficient intact");
  }

  { // Infinite max activity: x is unbounded above, so the row has no slack and
    // the whole row must be skipped.
    //
    // The rhs of -5 is what makes this fixture bite. The tempting refactor is
    // to *skip the unbounded term* and carry on with the rest of the row rather
    // than abandoning the row -- and that is unsafe: dropping x's contribution
    // gives a max activity of 0 and hence a slack of 5, shrinking -100 to -5
    // and turning the row into x - 5 y <= -5. That admits x = 0, y = 1 but no
    // longer x = 95, y = 1, which the original row does. So the assertion below
    // is on the coefficient surviving untouched, not merely on run() returning
    // false.
    Model mo;
    const int x = mo.addCol(0.0, COIN_DBL_MAX, false);
    const int y = mo.addCol(0.0, 1.0, true);
    mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, -5.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "infinite max activity: no change");
    checkClose(rowCoef(si, 0, y), -100.0, 1e-9, "infinite activity: coefficient intact");
  }

  { // singleton row: 100 y <= 50, slack 50 -- but only one coefficient
    Model mo;
    const int y = mo.addCol(0.0, 1.0, true);
    mo.addRow({ y }, { 100.0 }, -COIN_DBL_MAX, 50.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "singleton row: no change");
    checkClose(rowCoef(si, 0, y), 100.0, 1e-9, "singleton row: coefficient intact");
  }

  { // continuous oversized column: t is NOT integer, so it has no unit step
    Model mo;
    const int x = mo.addCol(0.0, 7.0, false);
    const int t = mo.addCol(0.0, 1.0, false);
    mo.addRow({ x, t }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "continuous column: no change");
    checkClose(rowCoef(si, 0, t), -100.0, 1e-9, "continuous column: coefficient intact");
  }

  { // fixed integer column: y in [1,1] has no step away from its bound. The
    // slack here is genuinely positive -- max activity is 7 - 100 = -93 against
    // an rhs of -100, so slack 7 -- which is what makes this a real test of the
    // fixed-column guard rather than of the slack guard.
    Model mo;
    const int x = mo.addCol(0.0, 7.0, false);
    const int y = mo.addCol(1.0, 1.0, true);
    mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, -100.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "fixed column: no change");
    checkClose(rowCoef(si, 0, y), -100.0, 1e-9, "fixed column: coefficient intact");
  }

  { // fractionally-bounded integer column: [0, 0.5] admits only y = 0, so there
    // is no unit step and the rhs adjustment would shave a sliver off the
    // feasible set. papilo asserts ub - lb >= 1 for exactly this reason.
    Model mo;
    const int x = mo.addCol(0.0, 7.0, false);
    const int y = mo.addCol(0.0, 0.5, true);
    mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0),
      "integer column with a sub-unit range [0,0.5]: no change");
    checkClose(rowCoef(si, 0, y), -100.0, 1e-9, "sub-unit range: coefficient intact");
  }

  { // non-positive slack: the row cannot be violated at all, so there is
    // nothing to strengthen even though 100 dwarfs the other coefficient
    Model mo;
    const int x = mo.addCol(0.0, 2.0, true);
    const int y = mo.addCol(0.0, 1.0, true);
    mo.addRow({ x, y }, { 1.0, 100.0 }, -COIN_DBL_MAX, 200.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "non-positive slack: no change");
    checkClose(rowCoef(si, 0, y), 100.0, 1e-9, "non-positive slack: coefficient intact");
  }
}

/* ------------------------------------------------- 6. non-integral slacks */

/* A slack that is genuinely fractional is used as-is; one that is a hair below
 * an integer is snapped UP to it. Rounding must only ever go up: a larger cap
 * is a weaker (safe) tightening, while a smaller one can cut off integer
 * points. This is the guard for the 105 non-integral-slack coefficients
 * measured across the mip-sanity corpus.
 *
 * Fixture: 2.5 x - 40 y <= 5,  x integer in [0,3],  y in {0,1}.
 * Max activity 7.5 (x at 3, y at 0), slack = 2.5 -- not integral and not within
 * tolerance of one, so -40 becomes exactly -2.5 and the rhs moves off y's lower
 * bound 0, i.e. not at all. Rounding the slack DOWN to 2 would be the bug, so
 * the assertion is on the exact value.
 */
void testNonIntegralSlack()
{
  printf("\n--- 6. genuinely fractional slack is used as-is ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 3.0, true, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0);
  mo.addRow({ x, y }, { 2.5, -40.0 }, -COIN_DBL_MAX, 5.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  const std::vector< std::string > feasibleBefore = enumerateFeasible(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);

  checkEqInt(cs.nCoefficients(), 1, "one coefficient tightened");
  checkClose(rowCoef(si, 0, y), -2.5, 1e-9,
    "coefficient of y is exactly -2.5 (slack NOT rounded down to 2)");
  check(enumerateFeasible(si) == feasibleBefore,
    "integer-feasible set unchanged (full enumeration)");
  check(!feasibleBefore.empty(), "enumeration found feasible points at all");
}

/* A slack that is a hair *below* an integer must snap up to that integer rather
 * than shave a sliver off the feasible set.
 *
 *     3 x - 40 y <= 6.0000001,   x integer in [0,3],  y in {0,1}
 *
 * Max activity 9, slack = 9 - 6.0000001 = 2.9999999, which is within the
 * feasibility tolerance of 3 and so becomes exactly 3. Leaving it at 2.9999999
 * would be a slightly-too-strong cap.
 */
void testSlackSnapsUp()
{
  printf("\n--- 6b. near-integral slack snaps UP to the integer ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 3.0, true, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0);
  mo.addRow({ x, y }, { 3.0, -40.0 }, -COIN_DBL_MAX, 6.0000001);

  OsiClpSolverInterface si;
  mo.loadInto(si);
  const std::vector< std::string > feasibleBefore = enumerateFeasible(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);

  checkEqInt(cs.nCoefficients(), 1, "one coefficient tightened");
  checkClose(rowCoef(si, 0, y), -3.0, 1e-9,
    "slack 2.9999999 snapped UP to exactly 3");
  check(rowCoef(si, 0, y) <= -3.0 + 1e-12,
    "coefficient magnitude not SMALLER than the true slack (would cut points)");
  check(enumerateFeasible(si) == feasibleBefore,
    "integer-feasible set unchanged (full enumeration)");
}

/* ----------------------------------------- 7. matrix views agree afterwards */

/* modifyCoefficient() does not invalidate the solver's cached row copy
 * (ClpModel::modifyCoefficient only clears whatsChanged_, never calling
 * freeCachedResults), so an implementation that pokes entries one at a time
 * leaves getMatrixByRow() reporting the OLD value while getMatrixByCol()
 * reports the new one -- and the LP then solves against a matrix nobody
 * intended. Batching into a copy and calling replaceMatrixOptional() once frees
 * the caches. This check fails loudly on the former.
 */
void testMatrixViewsAgree()
{
  printf("\n--- 7. row and column matrix views agree afterwards ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 7.0, false, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0);
  const int z = mo.addCol(0.0, 1.0, true, 1.0);
  mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
  mo.addRow({ x, z }, { 1.0, -60.0 }, -COIN_DBL_MAX, 0.0);
  mo.addRow({ x }, { 1.0 }, 7.0, 7.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);
  checkEqInt(cs.nCoefficients(), 2, "both big-M rows tightened");
  checkEqInt(cs.nRows(), 2, "two rows reported");
  checkClose(rowCoef(si, 0, y), -7.0, 1e-9, "row 0 reads the NEW value via the row view");
  checkClose(rowCoef(si, 1, z), -7.0, 1e-9, "row 1 reads the NEW value via the row view");

  CoinPackedMatrix fromCol;
  fromCol.reverseOrderedCopyOf(*si.getMatrixByCol());
  const CoinPackedMatrix *byRow = si.getMatrixByRow();
  checkEqInt(static_cast< int >(byRow->getNumElements()),
    static_cast< int >(fromCol.getNumElements()),
    "both views report the same number of elements");

  int mismatches = 0;
  for (int r = 0; r < si.getNumRows(); r++)
    for (int c = 0; c < si.getNumCols(); c++) {
      double a = 0.0, b = 0.0;
      {
        const CoinBigIndex s = byRow->getVectorStarts()[r];
        const int len = byRow->getVectorLengths()[r];
        for (CoinBigIndex k = s; k < s + len; k++)
          if (byRow->getIndices()[k] == c)
            a = byRow->getElements()[k];
      }
      {
        const CoinBigIndex s = fromCol.getVectorStarts()[r];
        const int len = fromCol.getVectorLengths()[r];
        for (CoinBigIndex k = s; k < s + len; k++)
          if (fromCol.getIndices()[k] == c)
            b = fromCol.getElements()[k];
      }
      if (fabs(a - b) > 1e-9)
        ++mismatches;
    }
  checkEqInt(mismatches, 0, "no row/column view disagreement (caches were freed)");

  // And the LP the solver actually solves must be the tightened one: y and z
  // are both forced to 1 now, so the objective is 2 rather than 0.07 + 7/60.
  si.initialSolve();
  check(si.isProvenOptimal(), "LP still solves to proven optimality afterwards");
  checkClose(si.getObjValue(), 2.0, 1e-6,
    "the solved LP is the TIGHTENED one (y = z = 1)");
}

/* ---------------------------------------------- 8. idempotence, one pass */

/* The slack is invariant under the change -- the coefficient loses exactly what
 * the rhs loses -- so a second pass must find nothing. That invariance is why
 * the production code runs a single pass with no iteration loop, and it was
 * confirmed on real instances (p0548/lseu/p0201/p2756 gain nothing from a
 * second pass).
 */
void testIdempotent()
{
  printf("\n--- 8. a second pass finds nothing (single pass suffices) ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 7.0, false, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0);
  mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
  mo.addRow({ x }, { 1.0 }, 7.0, 7.0);

  OsiClpSolverInterface si;
  mo.loadInto(si);

  CbcCoefficientStrengthening first;
  check(first.run(&si, si.messageHandler(), 0), "first pass changes something");
  const double afterFirst = lpBound(si);

  CbcCoefficientStrengthening second;
  check(!second.run(&si, si.messageHandler(), 0), "second pass changes nothing");
  checkEqInt(second.nCoefficients(), 0, "second pass: 0 coefficients");
  checkClose(lpBound(si), afterFirst, 1e-9, "LP bound unchanged by the second pass");
}

/* --------------------------------------------- 9. multi-row, mixed rows */

/* A model with one tightenable row among rows that must be left alone. The
 * aggregate counts alone would not catch a pass that fired on the wrong rows,
 * so every row's coefficient is asserted individually.
 */
void testMultiRow()
{
  printf("\n--- 9. multiple rows, mixed eligibility ---\n");
  Model mo;
  const int x = mo.addCol(0.0, 5.0, false, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0);
  const int z = mo.addCol(0.0, 1.0, true, 1.0);
  const int t = mo.addCol(0.0, 1.0, false, 0.0);
  mo.addRow({ x, y }, { 1.0, -80.0 }, -COIN_DBL_MAX, 0.0); // tightens: -80 -> -5
  mo.addRow({ x, z }, { 1.0, -80.0 }, -1.0, 0.0);          // range row: skipped
  mo.addRow({ x, t }, { 1.0, -80.0 }, -COIN_DBL_MAX, 0.0); // continuous: skipped
  mo.addRow({ y, z }, { 1.0, 1.0 }, -COIN_DBL_MAX, 1.0);   // nothing oversized

  OsiClpSolverInterface si;
  mo.loadInto(si);

  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);

  checkEqInt(cs.nCoefficients(), 1, "exactly one coefficient tightened");
  checkEqInt(cs.nRows(), 1, "exactly one row touched");
  checkClose(rowCoef(si, 0, y), -5.0, 1e-9, "row 0 (eligible): -80 -> -5");
  checkClose(rowCoef(si, 1, z), -80.0, 1e-9, "row 1 (range): intact");
  checkClose(rowCoef(si, 2, t), -80.0, 1e-9, "row 2 (continuous): intact");
  checkClose(rowCoef(si, 3, y), 1.0, 1e-9, "row 3 (nothing oversized): intact");
  checkClose(si.getRowLower()[1], -1.0, 1e-9, "row 1 lower bound intact");
  checkClose(si.getRowUpper()[1], 0.0, 1e-9, "row 1 upper bound intact");
}

/* ------------------------------------------------- 10. degenerate inputs */

void testDegenerate()
{
  printf("\n--- 10. degenerate inputs ---\n");
  {
    CbcCoefficientStrengthening cs;
    check(!cs.run(NULL, NULL, 0), "null solver: no crash, no change");
  }
  {
    OsiClpSolverInterface si;
    si.messageHandler()->setLogLevel(0);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "empty model: no crash, no change");
    checkEqInt(cs.nCoefficients(), 0, "empty model: 0 coefficients");
  }
  {
    Model mo; // a pure LP: no integer column anywhere
    const int a = mo.addCol(0.0, 7.0, false);
    const int b = mo.addCol(0.0, 1.0, false);
    mo.addRow({ a, b }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
    OsiClpSolverInterface si;
    mo.loadInto(si);
    CbcCoefficientStrengthening cs;
    check(!cs.run(&si, si.messageHandler(), 0), "pure LP: no change");
  }
}

/* ------------------------------- 11. end to end through -coefStrengthening */

/* Everything above drives the class directly. This drives the whole solver, so
 * it also covers the wiring: the parameter, the call site in
 * preRootLPStrenghtening(), and the fact that the pass runs before the root LP.
 *
 * The observable outcome must be identical either way -- the optimum is 1
 * regardless -- and that is itself the assertion worth making: a strengthening
 * step that changed the answer would be a bug, not a feature. Preprocessing is
 * switched off so nothing else tightens the row on the way through. The bound
 * difference the pass produces is asserted at the class level below, since the
 * root LP value is not exposed through the C interface.
 */
void testParameterEndToEnd()
{
  printf("\n--- 11. end to end: -coefStrengthening on/off ---\n");

  for (int off = 0; off < 2; off++) {
    Cbc_Model *model = Cbc_newModel();
    Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 0);
    if (off)
      Cbc_setParameter(model, "coefStrengthening", "off");
    Cbc_setParameter(model, "preprocess", "off");

    Cbc_addCol(model, "x", 0.0, 7.0, 0.0, 0, 0, NULL, NULL);
    Cbc_addCol(model, "y", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    { // x - 100 y <= 0
      int cols[2] = { 0, 1 };
      double coefs[2] = { 1.0, -100.0 };
      Cbc_addRow(model, "bigM", 2, cols, coefs, 'L', 0.0);
    }
    { // x = 7
      int cols[1] = { 0 };
      double coefs[1] = { 1.0 };
      Cbc_addRow(model, "pin", 1, cols, coefs, 'E', 7.0);
    }

    Cbc_solve(model);
    checkClose(Cbc_getObjValue(model), 1.0, 1e-6,
      off ? "coefStrengthening off: optimum is 1 (y must be 1)"
          : "coefStrengthening on: optimum is 1 (y must be 1) -- unchanged answer");
    Cbc_deleteModel(model);
  }

  // The root LP the pass hands on is the part that differs, asserted directly.
  Model mo;
  const int x = mo.addCol(0.0, 7.0, false, 0.0);
  const int y = mo.addCol(0.0, 1.0, true, 1.0);
  mo.addRow({ x, y }, { 1.0, -100.0 }, -COIN_DBL_MAX, 0.0);
  mo.addRow({ x }, { 1.0 }, 7.0, 7.0);
  OsiClpSolverInterface si;
  mo.loadInto(si);
  const double lpOff = lpBound(si);
  CbcCoefficientStrengthening cs;
  cs.run(&si, si.messageHandler(), 0);
  const double lpOn = lpBound(si);
  check(lpOn > lpOff,
    "the root LP the pass hands on is strictly tighter than the untouched one");
  printf("        (root LP bound: untouched = %.6g, tightened = %.6g)\n", lpOff, lpOn);
}

} // end anonymous namespace

int main(int argc, char **argv)
{
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      printf("coefstr-test — tests for CbcCoefficientStrengthening (coefficient\n"
             "tightening / \"big-M strengthening\", step 3 of the pre-root-LP\n"
             "strengthening phase).\n\n"
             "Usage: coefstr-test\n\n"
             "Takes no options and needs no data files: every fixture and its\n"
             "expected coefficients and right-hand sides are built in-process and\n"
             "worked out by hand in the comment above each test.\n\n"
             "Exit code: 0 = all checks passed, 1 = one or more failed.\n");
      return 0;
    }
  }

  printf("=== coefstr-test: coefficient tightening (big-M strengthening) ===\n");

  testTextbookBigM();
  testGeneralInteger();
  testNegativeCoefficient();
  testGreaterEqualRow();
  testMustNotFire();
  testNonIntegralSlack();
  testSlackSnapsUp();
  testMatrixViewsAgree();
  testIdempotent();
  testMultiRow();
  testDegenerate();
  testParameterEndToEnd();

  printf("\n");
  if (g_failures) {
    printf("RESULT: %d check(s) FAILED\n", g_failures);
    return 1;
  }
  printf("RESULT: all checks passed\n");
  return 0;
}
