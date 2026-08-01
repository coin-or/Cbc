/*
 * lazy-constraint-test — targeted tests for the C interface's lazy-constraint
 * machinery (Cbc_addLazyConstraint, and Cbc_addCutCallback with
 * atSolution=1), driven by the TSP with subtour elimination. That is the
 * canonical use of the feature: the SEC family has exponentially many members,
 * so it is stated lazily and only the violated ones ever reach a solver.
 *
 * Why the TSP rather than a two-variable toy. A lazy constraint that is
 * silently *ignored* still yields an optimal-looking answer -- just the answer
 * to a relaxation. The assignment problem (degree constraints only) is exactly
 * the TSP minus its SECs, and it is a strictly weaker model on every instance
 * used here, so the two are told apart by the objective alone. The reference
 * optimum is not hardcoded: each instance is built by a fixed LCG and solved
 * three independent ways -- Held-Karp dynamic programming in this file, the
 * same MIP with SECs as *regular rows*, and the same MIP with SECs as lazy
 * constraints -- which must agree, while the SEC-free model must come out
 * strictly below. Sizes stay at n <= 12 so full SEC enumeration and exact DP
 * are both cheap (the whole program runs in a few seconds).
 *
 * Coverage:
 *
 *   1. Enforcement. Lazy SECs must bind exactly as row SECs do: same optimum
 *      as Held-Karp on every instance, and strictly above the SEC-free bound.
 *      The *published* solution is also walked as a successor array, so it must
 *      be integral and one Hamiltonian cycle of exactly that cost -- an
 *      objective check alone would pass on a solution that is two subtours
 *      priced as if it were a tour.
 *   2. LP relaxation untouched. A lazy constraint is not part of the model the
 *      LP solves, so Cbc_solveLinearProgram() must return the same bound with
 *      and without them, and Cbc_resolve() must not drift from it.
 *   3. Repeated Cbc_solve() on one model. Each solve must independently reach
 *      the optimum. Cbc reuses the previous solve's answer as a MIPStart, and
 *      a MIPStart is screened *before* branchAndBound() marks the model as
 *      carrying lazy constraints -- so this is where a screen that reads the
 *      wrong solution vector shows up, as a model with a perfectly good
 *      optimum reporting proven infeasibility (objective 1.8e308).
 *   4. Lazy constraints added *after* a first solve. The second solve must move
 *      from the SEC-free bound to the true optimum: constraints registered
 *      between solves have to reach the search, and a first-solve answer that
 *      violates them must not survive as an incumbent.
 *   5. Cbc_clone(). Lazy constraints are part of the model's definition, not a
 *      cache of a previous solve, so a clone must carry them -- a clone that
 *      dropped them would solve the assignment relaxation and report it as
 *      proven optimal. Also checked to be a deep copy: a constraint added to
 *      the clone must not bind the original.
 *   6. MIPStart. An externally supplied solution that *violates* a lazy
 *      constraint must be rejected rather than trusted: taking it as an
 *      incumbent installs a cutoff that can exclude the genuine optimum.
 *      Checked in both senses, since the cutoff is negated for maximization.
 *      A merely suboptimal but lazy-feasible MIPStart must be accepted and
 *      improved upon.
 *   7. Cbc_addCutCallback(..., atSolution=1) — the other way to express a lazy
 *      family, separating on demand from the integral solution instead of
 *      enumerating up front. Must reach the same optimum.
 *   8. The iterative LP cut loop (Cbc_solveLinearProgram / Cbc_addRow /
 *      Cbc_resolve): bound must increase monotonically toward the SEC-closure
 *      value, and the MIP solved afterwards must still reach the optimum.
 *   9. Staged rows and columns reaching the solve. Cbc_addCol()/Cbc_addRow()
 *      only append to staging arrays, so a solve that forgets to flush them
 *      reads an empty or stale solver: a fresh model gets dispatched as an LP
 *      because it appears to have no integers (returning the *relaxation* while
 *      reporting proven optimality), and a row added after an earlier solve is
 *      dropped from the model it was meant to constrain.
 *  10. Solution pool. Lazy constraints promise that *no* accepted solution
 *      violates one, so every solution in the pool -- not just the best -- must
 *      be feasible and a single Hamiltonian cycle.
 *  11. A solution offered from outside the search -- an explicit MIPStart or
 *      the previous answer replayed on a re-solve -- both screened against the
 *      right vector and confirmed at its own bounds. Screening means re-running
 *      the at-solution generators, which read the *solver's* column solution
 *      while the candidate lives elsewhere, so the candidate has to be installed
 *      first; and the winner is confirmed by fixing its integer variables and
 *      re-optimising the continuous ones, so that fixing has to survive to the
 *      resolve. Reading the wrong vector rejects the optimum and reports proven
 *      infeasibility; losing the fixing publishes the root LP relaxation
 *      instead, at a point that is fractional or (after rounding) violates the
 *      model's own row. The TSP cases above cannot see either: the vector left
 *      in the solver is the previous optimal tour, so it gives the right verdict
 *      by accident, and the SEC-closed relaxation is tight at the optimum. Two
 *      binaries under x + y <= 1.5 leave a fractional vector (x = y = 0.75) and
 *      a strictly loose relaxation, exposing both. The published point is
 *      checked, not just the objective.
 *
 * Usage: lazy-constraint-test [--help]
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "Cbc_C_Interface.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
    printf("  ok:   %s (actual=%.10g expected=%.10g delta=%.3e <= tol=%.3e)\n",
      msg.c_str(), actual, expected, delta, tol);
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

/* ── TSP instance ─────────────────────────────────────────────────────── */

/* A symmetric distance matrix from a fixed linear congruential generator, so
 * the instance is byte-identical on every platform and no fixture file is
 * needed. Distances land in [10, 259], comfortably away from degeneracy. */
struct Tsp {
  int n;
  std::vector<int> d;

  Tsp(int nn, unsigned seed)
    : n(nn), d(static_cast<size_t>(nn) * nn, 0)
  {
    unsigned s = seed;
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        s = s * 1103515245u + 12345u;
        int v = 10 + static_cast<int>((s >> 16) % 250);
        d[i * n + j] = v;
        d[j * n + i] = v;
      }
    }
  }

  int dist(int i, int j) const { return d[i * n + j]; }
  /* Arc variable x(i,j) is column i*n+j. The n diagonal columns exist but are
   * fixed to zero, which buys trivial index arithmetic for the price of n
   * columns the presolve removes anyway. */
  int col(int i, int j) const { return i * n + j; }
};

/* Exact TSP optimum by Held-Karp dynamic programming over subsets of
 * {1..n-1}. Independent of Cbc, and the reference every MIP result below is
 * checked against. */
double heldKarp(const Tsp &t)
{
  const int n = t.n;
  const int full = 1 << (n - 1);
  const double INF = 1e30;
  std::vector<double> dp(static_cast<size_t>(full) * n, INF);
  for (int j = 1; j < n; ++j)
    dp[static_cast<size_t>(1 << (j - 1)) * n + j] = t.dist(0, j);
  for (int mask = 0; mask < full; ++mask) {
    for (int j = 1; j < n; ++j) {
      double base = dp[static_cast<size_t>(mask) * n + j];
      if (!(mask & (1 << (j - 1))) || base >= INF)
        continue;
      for (int k = 1; k < n; ++k) {
        if (mask & (1 << (k - 1)))
          continue;
        int nmask = mask | (1 << (k - 1));
        double v = base + t.dist(j, k);
        if (v < dp[static_cast<size_t>(nmask) * n + k])
          dp[static_cast<size_t>(nmask) * n + k] = v;
      }
    }
  }
  double best = INF;
  for (int j = 1; j < n; ++j) {
    double v = dp[static_cast<size_t>(full - 1) * n + j] + t.dist(j, 0);
    if (v < best)
      best = v;
  }
  return best;
}

enum SecMode {
  SEC_NONE, /* assignment relaxation: degree constraints only */
  SEC_ROWS, /* SECs as ordinary rows -- the reference formulation */
  SEC_LAZY  /* SECs via Cbc_addLazyConstraint -- what is under test */
};

/* Directed assignment formulation plus, per mode, the full SEC family for
 * every 2 <= |S| <= n-1. Exponential by construction, which is the point: it
 * is only tractable here because n <= 12.
 * `sense` is +1 to minimize the tour cost or -1 to maximize its negation --
 * the same problem, reached through Cbc's other objective sense. */
Cbc_Model *buildTsp(const Tsp &t, SecMode mode, double sense = 1.0)
{
  const int n = t.n;
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, sense);
  Cbc_setLogLevel(m, 0);

  char nm[32];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      snprintf(nm, sizeof(nm), "x(%d,%d)", i, j);
      Cbc_addCol(m, nm, 0.0, (i == j) ? 0.0 : 1.0,
        sense * t.dist(i, j), 1, 0, NULL, NULL);
    }
  }

  std::vector<int> idx(static_cast<size_t>(n) * n);
  std::vector<double> cf(static_cast<size_t>(n) * n, 1.0);

  for (int i = 0; i < n; ++i) {
    int k = 0;
    for (int j = 0; j < n; ++j)
      if (j != i)
        idx[k++] = t.col(i, j);
    snprintf(nm, sizeof(nm), "out(%d)", i);
    Cbc_addRow(m, nm, k, &idx[0], &cf[0], 'E', 1.0);
  }
  for (int j = 0; j < n; ++j) {
    int k = 0;
    for (int i = 0; i < n; ++i)
      if (i != j)
        idx[k++] = t.col(i, j);
    snprintf(nm, sizeof(nm), "in(%d)", j);
    Cbc_addRow(m, nm, k, &idx[0], &cf[0], 'E', 1.0);
  }

  if (mode != SEC_NONE) {
    for (int mask = 1; mask < (1 << n); ++mask) {
      int sz = 0;
      for (int b = 0; b < n; ++b)
        if (mask & (1 << b))
          ++sz;
      if (sz < 2 || sz > n - 1)
        continue;
      int k = 0;
      for (int i = 0; i < n; ++i) {
        if (!(mask & (1 << i)))
          continue;
        for (int j = 0; j < n; ++j)
          if (i != j && (mask & (1 << j)))
            idx[k++] = t.col(i, j);
      }
      if (mode == SEC_LAZY)
        Cbc_addLazyConstraint(m, k, &idx[0], &cf[0], 'L', sz - 1);
      else
        Cbc_addRow(m, "sec", k, &idx[0], &cf[0], 'L', sz - 1);
    }
  }
  return m;
}

/* Walk the solution as a successor array from node 0. Returns true only for an
 * integral vector that is a single Hamiltonian cycle; `cost` gets the tour cost
 * recomputed from the distance matrix, independently of the objective Cbc
 * reports, and `nFrac` the count of fractional arcs. Two disjoint subtours
 * satisfy every degree constraint and would pass an objective-only check, so
 * this is what actually distinguishes an enforced SEC family from an ignored
 * one at the level of the published solution. */
bool isSingleTour(const Tsp &t, const double *x, double &cost, int &nFrac)
{
  const int n = t.n;
  std::vector<int> succ(n, -1);
  cost = 0.0;
  nFrac = 0;
  int outDeg = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double v = x[t.col(i, j)];
      if (v > 1.0e-6 && v < 1.0 - 1.0e-6)
        ++nFrac;
      if (v > 0.5) {
        if (succ[i] >= 0)
          return false; /* out-degree above one: not a successor function */
        succ[i] = j;
        cost += t.dist(i, j);
        ++outDeg;
      }
    }
  }
  if (nFrac != 0 || outDeg != n)
    return false;
  int at = 0, len = 0;
  do {
    at = succ[at];
    ++len;
  } while (at != 0 && len <= n);
  return at == 0 && len == n;
}

/* One combined verdict on a solved model: objective as reported, and the
 * published solution walked as a tour of exactly that cost. */
void checkTourSolution(Cbc_Model *m, const Tsp &t, double expectedObj,
  const std::string &label)
{
  checkClose(Cbc_getObjValue(m), expectedObj, 1.0e-6, label + ": objective");
  const double *x = Cbc_getColSolution(m);
  if (!x) {
    ++g_failures;
    printf("  FAIL: %s: no solution published\n", label.c_str());
    return;
  }
  double cost = 0.0;
  int nFrac = 0;
  bool tour = isSingleTour(t, x, cost, nFrac);
  check(tour, label + ": published solution is one Hamiltonian cycle (frac="
    + std::to_string(nFrac) + ")");
  checkClose(cost, fabs(expectedObj), 1.0e-6,
    label + ": tour cost recomputed from the distance matrix");
}

/* ── Test 1: lazy SECs bind exactly as row SECs do ─────────────────────── */

void testEnforcement(const Tsp &t)
{
  printf("\n=== Test 1: lazy SECs enforce the tour constraint — n=%d ===\n", t.n);
  double hk = heldKarp(t);
  printf("  Held-Karp exact optimum = %.0f\n", hk);

  Cbc_Model *none = buildTsp(t, SEC_NONE);
  Cbc_solve(none);
  double assignBound = Cbc_getObjValue(none);
  check(Cbc_isProvenOptimal(none) != 0, "assignment relaxation solved to optimality");
  /* The whole test rests on this: if the SEC-free model already priced at the
   * tour optimum, an ignored SEC family would be indistinguishable. */
  check(assignBound < hk - 1.0e-6, "assignment bound " + std::to_string(assignBound)
    + " is strictly below the tour optimum (so SECs are discriminating)");
  Cbc_deleteModel(none);

  Cbc_Model *rows = buildTsp(t, SEC_ROWS);
  Cbc_solve(rows);
  check(Cbc_isProvenOptimal(rows) != 0, "SECs as rows: proven optimal");
  checkTourSolution(rows, t, hk, "SECs as rows");
  Cbc_deleteModel(rows);

  Cbc_Model *lazy = buildTsp(t, SEC_LAZY);
  Cbc_solve(lazy);
  check(Cbc_isProvenOptimal(lazy) != 0, "SECs as lazy constraints: proven optimal");
  check(Cbc_isProvenInfeasible(lazy) == 0,
    "SECs as lazy constraints: not reported infeasible");
  checkTourSolution(lazy, t, hk, "SECs as lazy constraints");
  Cbc_deleteModel(lazy);
}

/* ── Test 2: the LP relaxation must not see them ───────────────────────── */

void testLpRelaxationUnaffected(const Tsp &t)
{
  printf("\n=== Test 2: lazy constraints do not tighten the LP relaxation — n=%d ===\n",
    t.n);

  Cbc_Model *plain = buildTsp(t, SEC_NONE);
  Cbc_solveLinearProgram(plain);
  double lpPlain = Cbc_getObjValue(plain);
  Cbc_deleteModel(plain);

  Cbc_Model *lazy = buildTsp(t, SEC_LAZY);
  Cbc_solveLinearProgram(lazy);
  checkClose(Cbc_getObjValue(lazy), lpPlain, 1.0e-6,
    "LP relaxation identical with and without lazy constraints");
  Cbc_resolve(lazy);
  checkClose(Cbc_getObjValue(lazy), lpPlain, 1.0e-6,
    "Cbc_resolve() does not drift from it either");
  Cbc_deleteModel(lazy);
}

/* ── Test 3: repeated Cbc_solve() on one model ─────────────────────────── */

void testRepeatedSolve(const Tsp &t)
{
  printf("\n=== Test 3: three consecutive Cbc_solve() on one lazy model — n=%d ===\n",
    t.n);
  /* Solves after the first inherit the previous answer as a MIPStart, which is
   * applied and screened before branchAndBound() flags the model as carrying
   * lazy constraints. A screen reading the wrong solution vector in that window
   * rejects a perfectly good incumbent and, because the objective is bound by
   * reference, corrupts it into "proven infeasible" (1.8e308). */
  double hk = heldKarp(t);
  Cbc_Model *m = buildTsp(t, SEC_LAZY);
  for (int s = 1; s <= 3; ++s) {
    Cbc_solve(m);
    std::string label = "solve " + std::to_string(s);
    check(Cbc_isProvenOptimal(m) != 0, label + ": proven optimal");
    check(Cbc_isProvenInfeasible(m) == 0, label + ": not reported infeasible");
    checkTourSolution(m, t, hk, label);
  }
  Cbc_deleteModel(m);
}

/* ── Test 4: lazy constraints registered between solves ────────────────── */

void testAddAfterSolve(const Tsp &t)
{
  printf("\n=== Test 4: lazy constraints added after a first solve — n=%d ===\n", t.n);
  double hk = heldKarp(t);

  Cbc_Model *m = buildTsp(t, SEC_NONE);
  Cbc_solve(m);
  double first = Cbc_getObjValue(m);
  check(first < hk - 1.0e-6,
    "first solve (no SECs) returns the weaker assignment bound "
      + std::to_string(first));

  /* Now state the SEC family and solve again. The first solve's answer, which
   * violates it, is carried in as a MIPStart; accepting it on trust would
   * install a cutoff below the true optimum and the search would report the
   * model infeasible. */
  const int n = t.n;
  std::vector<int> idx(static_cast<size_t>(n) * n);
  std::vector<double> cf(static_cast<size_t>(n) * n, 1.0);
  int nSec = 0;
  for (int mask = 1; mask < (1 << n); ++mask) {
    int sz = 0;
    for (int b = 0; b < n; ++b)
      if (mask & (1 << b))
        ++sz;
    if (sz < 2 || sz > n - 1)
      continue;
    int k = 0;
    for (int i = 0; i < n; ++i) {
      if (!(mask & (1 << i)))
        continue;
      for (int j = 0; j < n; ++j)
        if (i != j && (mask & (1 << j)))
          idx[k++] = t.col(i, j);
    }
    Cbc_addLazyConstraint(m, k, &idx[0], &cf[0], 'L', sz - 1);
    ++nSec;
  }
  printf("  %d lazy SECs registered after the first solve\n", nSec);

  Cbc_solve(m);
  check(Cbc_isProvenOptimal(m) != 0, "second solve: proven optimal");
  check(Cbc_isProvenInfeasible(m) == 0, "second solve: not reported infeasible");
  checkTourSolution(m, t, hk, "second solve");
  Cbc_deleteModel(m);
}

/* ── Test 5: Cbc_clone() carries them, and deeply ──────────────────────── */

void testClone(const Tsp &t)
{
  printf("\n=== Test 5: Cbc_clone() and lazy constraints — n=%d ===\n", t.n);
  double hk = heldKarp(t);

  Cbc_Model *orig = buildTsp(t, SEC_LAZY);
  Cbc_Model *cl = Cbc_clone(orig);
  Cbc_setLogLevel(cl, 0);

  /* Forbid one arbitrary Hamiltonian cycle -- 0->1->...->n-1->0, in both
   * directions -- in the clone only. The clone's optimum may then be no better
   * than the original's, and the original must be unmoved. */
  const int n = t.n;
  std::vector<int> idx;
  std::vector<double> cf;
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    idx.push_back(t.col(i, j));
    cf.push_back(1.0);
    idx.push_back(t.col(j, i));
    cf.push_back(1.0);
  }
  Cbc_addLazyConstraint(cl, static_cast<int>(idx.size()), &idx[0], &cf[0], 'L',
    n - 1);

  Cbc_solve(cl);
  check(Cbc_isProvenOptimal(cl) != 0, "clone: proven optimal");
  /* The point of the check: a clone that dropped the SECs solves the
   * assignment relaxation and reports it as optimal, i.e. comes out *below*
   * the tour optimum. */
  double clObj = Cbc_getObjValue(cl);
  check(clObj >= hk - 1.0e-6, "clone kept the lazy constraints (obj "
    + std::to_string(clObj) + " >= tour optimum " + std::to_string(hk) + ")");
  double cost = 0.0;
  int nFrac = 0;
  check(isSingleTour(t, Cbc_getColSolution(cl), cost, nFrac),
    "clone: published solution is one Hamiltonian cycle");
  Cbc_deleteModel(cl);

  Cbc_solve(orig);
  check(Cbc_isProvenOptimal(orig) != 0, "original after cloning: proven optimal");
  checkTourSolution(orig, t, hk,
    "original unaffected by the constraint added to the clone");
  Cbc_deleteModel(orig);
}

/* ── Test 6: MIPStart screened against the lazy constraints ────────────── */

/* Two disjoint cycles over the first and second half of the nodes: every
 * degree constraint is satisfied, so nothing but a SEC rules it out. */
void setSubtourMipStart(Cbc_Model *m, const Tsp &t, double &startCost)
{
  const int n = t.n;
  const int h = n / 2;
  std::vector<std::string> names;
  std::vector<double> vals;
  startCost = 0.0;
  for (int i = 0; i < h; ++i) {
    int j = (i + 1) % h;
    char nm[32];
    snprintf(nm, sizeof(nm), "x(%d,%d)", i, j);
    names.push_back(nm);
    vals.push_back(1.0);
    startCost += t.dist(i, j);
  }
  for (int i = 0; i < n - h; ++i) {
    int j = (i + 1) % (n - h);
    char nm[32];
    snprintf(nm, sizeof(nm), "x(%d,%d)", h + i, h + j);
    names.push_back(nm);
    vals.push_back(1.0);
    startCost += t.dist(h + i, h + j);
  }
  std::vector<const char *> ptrs;
  for (size_t k = 0; k < names.size(); ++k)
    ptrs.push_back(names[k].c_str());
  Cbc_setMIPStart(m, static_cast<int>(ptrs.size()), &ptrs[0], &vals[0]);
}

/* The identity cycle 0->1->...->n-1->0: lazy-feasible, generally suboptimal. */
void setTourMipStart(Cbc_Model *m, const Tsp &t, double &startCost)
{
  const int n = t.n;
  std::vector<std::string> names;
  std::vector<double> vals;
  startCost = 0.0;
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    char nm[32];
    snprintf(nm, sizeof(nm), "x(%d,%d)", i, j);
    names.push_back(nm);
    vals.push_back(1.0);
    startCost += t.dist(i, j);
  }
  std::vector<const char *> ptrs;
  for (size_t k = 0; k < names.size(); ++k)
    ptrs.push_back(names[k].c_str());
  Cbc_setMIPStart(m, static_cast<int>(ptrs.size()), &ptrs[0], &vals[0]);
}

void testMipStart(const Tsp &t)
{
  printf("\n=== Test 6: MIPStart screened against the lazy constraints — n=%d ===\n",
    t.n);
  double hk = heldKarp(t);

  /* Minimize. A MIPStart is applied by CbcSolver before branchAndBound() runs,
   * so it is the one solution the search never derives itself and the one place
   * a lazy constraint can be bypassed outright. */
  {
    Cbc_Model *m = buildTsp(t, SEC_LAZY);
    double startCost = 0.0;
    setSubtourMipStart(m, t, startCost);
    printf("  minimize: MIPStart is two disjoint cycles, cost %.0f "
           "(degree-feasible, SEC-violating)\n", startCost);
    Cbc_solve(m);
    check(Cbc_isProvenOptimal(m) != 0, "minimize + SEC-violating MIPStart: proven optimal");
    check(Cbc_isProvenInfeasible(m) == 0,
      "minimize + SEC-violating MIPStart: not reported infeasible");
    checkTourSolution(m, t, hk, "minimize + SEC-violating MIPStart");
    Cbc_deleteModel(m);
  }

  /* Maximize the negated cost: the same problem, but the incumbent's cutoff is
   * installed negated, which is a separate code path. */
  {
    Cbc_Model *m = buildTsp(t, SEC_LAZY, -1.0);
    double startCost = 0.0;
    setSubtourMipStart(m, t, startCost);
    Cbc_solve(m);
    check(Cbc_isProvenOptimal(m) != 0, "maximize + SEC-violating MIPStart: proven optimal");
    checkTourSolution(m, t, -hk, "maximize + SEC-violating MIPStart");
    Cbc_deleteModel(m);
  }

  /* A lazy-feasible but suboptimal start must be kept and improved on, not
   * rejected: the screens above must not fire on a solution that satisfies
   * every lazy constraint. */
  {
    Cbc_Model *m = buildTsp(t, SEC_LAZY);
    double startCost = 0.0;
    setTourMipStart(m, t, startCost);
    check(startCost > hk + 1.0e-6, "identity-cycle MIPStart cost "
      + std::to_string(startCost) + " is suboptimal (so improving on it is visible)");
    Cbc_solve(m);
    check(Cbc_isProvenOptimal(m) != 0, "lazy-feasible suboptimal MIPStart: proven optimal");
    checkTourSolution(m, t, hk, "lazy-feasible suboptimal MIPStart");
    Cbc_deleteModel(m);
  }
}

/* ── Test 7: the same family via a cut callback at integer solutions ───── */

/* Parse "x(i,j)". Names, not column indices, are the only stable handle inside
 * a callback: presolve renumbers columns. */
bool parseArcName(const char *name, int &u, int &v)
{
  const char *s = strstr(name, "x(");
  if (!s)
    return false;
  s += 2;
  char *end = NULL;
  u = static_cast<int>(strtol(s, &end, 10));
  if (*end != ',')
    return false;
  v = static_cast<int>(strtol(end + 1, &end, 10));
  return *end == ')';
}

/* Separator: on an integral solution, find the smallest connected component of
 * the arc support and, if it is not the whole graph, add its SEC as a global
 * row cut. Registered with atSolution=1, which is what makes it a lazy
 * constraint rather than a strengthening cut. */
void subtourCutCallback(void *osiSolver, void *osiCuts, void *appData,
  int /*level*/, int /*pass*/)
{
  if (!Osi_isProvenOptimal(osiSolver))
    return;
  const int n = *static_cast<int *>(appData);
  const int ncols = Osi_getNumCols(osiSolver);
  const double *x = Osi_getColSolution(osiSolver);

  /* Map (i,j) -> current column index, by name. */
  std::vector<int> col(static_cast<size_t>(n) * n, -1);
  int nFrac = 0;
  for (int c = 0; c < ncols; ++c) {
    char nm[64] = "";
    Osi_getColName(osiSolver, c, nm, sizeof(nm));
    int u = 0, v = 0;
    if (!parseArcName(nm, u, v))
      continue;
    if (u < 0 || u >= n || v < 0 || v >= n)
      continue;
    col[u * n + v] = c;
    if (x[c] > 1.0e-6 && x[c] < 1.0 - 1.0e-6)
      ++nFrac;
  }
  if (nFrac)
    return; /* only separate on integer solutions */

  std::vector<char> seen(n, 0);
  std::vector<int> best, comp, queue;
  size_t bestSize = static_cast<size_t>(n) + 1;
  for (int start = 0; start < n; ++start) {
    if (seen[start])
      continue;
    comp.clear();
    queue.clear();
    queue.push_back(start);
    seen[start] = 1;
    for (size_t qh = 0; qh < queue.size(); ++qh) {
      int v = queue[qh];
      comp.push_back(v);
      for (int w = 0; w < n; ++w) {
        if (w == v || seen[w])
          continue;
        int a = col[v * n + w], b = col[w * n + v];
        /* Undirected connectivity: either orientation joins the component. */
        if ((a >= 0 && x[a] > 0.5) || (b >= 0 && x[b] > 0.5)) {
          seen[w] = 1;
          queue.push_back(w);
        }
      }
    }
    if (comp.size() < bestSize) {
      bestSize = comp.size();
      best = comp;
    }
  }
  if (static_cast<int>(bestSize) >= n)
    return; /* one component: a genuine Hamiltonian cycle, nothing to cut */

  std::vector<char> inS(n, 0);
  for (size_t k = 0; k < best.size(); ++k)
    inS[best[k]] = 1;
  /* Every arc inside S must appear, including the ones currently at zero --
   * the cut has to forbid the whole family, not just the arcs in use. */
  std::vector<int> cutIdx;
  std::vector<double> cutCoef;
  for (int u = 0; u < n; ++u) {
    if (!inS[u])
      continue;
    for (int v = 0; v < n; ++v) {
      if (u == v || !inS[v] || col[u * n + v] < 0)
        continue;
      cutIdx.push_back(col[u * n + v]);
      cutCoef.push_back(1.0);
    }
  }
  if (!cutIdx.empty())
    OsiCuts_addGlobalRowCut(osiCuts, static_cast<int>(cutIdx.size()), &cutIdx[0],
      &cutCoef[0], 'L', static_cast<double>(bestSize) - 1.0);
}

void testCutCallbackAtSolution(const Tsp &t)
{
  printf("\n=== Test 7: SECs via Cbc_addCutCallback(atSolution=1) — n=%d ===\n", t.n);
  double hk = heldKarp(t);
  int n = t.n;
  Cbc_Model *m = buildTsp(t, SEC_NONE);
  Cbc_addCutCallback(m, subtourCutCallback, "TSP-SEC", &n, 1, 1);
  Cbc_solve(m);
  check(Cbc_isProvenOptimal(m) != 0, "cut callback at solution: proven optimal");
  checkTourSolution(m, t, hk, "cut callback at solution");
  Cbc_deleteModel(m);
}

/* ── Test 8: the iterative LP cut loop ────────────────────────────────── */

void testLpCutLoop(const Tsp &t)
{
  printf("\n=== Test 8: iterative LP cut loop (solveLinearProgram/addRow/resolve) — n=%d ===\n",
    t.n);
  const int n = t.n;
  const int maxIter = 100;
  double hk = heldKarp(t);
  Cbc_Model *m = buildTsp(t, SEC_NONE);

  double first = 0.0, prev = -1.0e30;
  int iter = 0, added = 0, nonMonotone = 0;
  bool connected = false;
  while (iter < maxIter) {
    int rc = (iter == 0) ? Cbc_solveLinearProgram(m) : Cbc_resolve(m);
    ++iter;
    if (rc != 0) {
      ++g_failures;
      printf("  FAIL: LP solve returned status %d at iteration %d\n", rc, iter);
      break;
    }
    double obj = Cbc_getObjValue(m);
    if (iter == 1)
      first = obj;
    else if (obj < prev - 1.0e-6)
      ++nonMonotone;
    prev = obj;

    /* Connected components of the *support* of the LP solution -- every arc
     * carrying anything at all, not just the ones above 0.5. Each component S
     * has zero flow crossing its boundary, so its SEC is violated by the full
     * degree count; thresholding at 0.5 instead would miss that on a fractional
     * solution and could isolate single nodes, whose SEC is vacuous. When the
     * support is connected there is no zero-crossing cut left to find and the
     * loop is done (going further would need a genuine min-cut separator). */
    const double *x = Cbc_getColSolution(m);
    std::vector<char> seen(n, 0);
    std::vector<std::vector<int> > comps;
    for (int start = 0; start < n; ++start) {
      if (seen[start])
        continue;
      std::vector<int> comp, queue;
      queue.push_back(start);
      seen[start] = 1;
      for (size_t qh = 0; qh < queue.size(); ++qh) {
        int v = queue[qh];
        comp.push_back(v);
        for (int w = 0; w < n; ++w)
          if (!seen[w] && (x[t.col(v, w)] > 1.0e-6 || x[t.col(w, v)] > 1.0e-6)) {
            seen[w] = 1;
            queue.push_back(w);
          }
      }
      comps.push_back(comp);
    }
    if (comps.size() == 1u) {
      connected = true;
      break;
    }
    for (size_t c = 0; c < comps.size(); ++c) {
      const std::vector<int> &comp = comps[c];
      if (comp.size() < 2u || comp.size() > static_cast<size_t>(n) - 1u)
        continue; /* only 2 <= |S| <= n-1 gives a valid SEC */
      std::vector<int> idx;
      std::vector<double> cf;
      for (size_t a = 0; a < comp.size(); ++a)
        for (size_t b = 0; b < comp.size(); ++b)
          if (a != b) {
            idx.push_back(t.col(comp[a], comp[b]));
            cf.push_back(1.0);
          }
      Cbc_addRow(m, "sec", static_cast<int>(idx.size()), &idx[0], &cf[0], 'L',
        static_cast<double>(comp.size()) - 1.0);
      ++added;
    }
  }
  printf("  %d LP iterations, %d cut(s) added, bound %.0f -> %.0f\n",
    iter, added, first, prev);
  check(nonMonotone == 0, "LP bound never decreased while adding valid cuts");
  check(connected, "cut loop ended with a connected support, not at the "
    + std::to_string(maxIter) + "-iteration cap");
  check(added > 0, "the loop actually separated something ("
    + std::to_string(added) + " cut(s))");
  check(prev > first + 1.0e-6, "the added rows reached the LP: bound rose from "
    + std::to_string(first) + " to " + std::to_string(prev));
  /* Every added row is a valid SEC, so the LP is still a relaxation of the
   * TSP: a bound above the true optimum would mean one of them cut it off. */
  check(prev <= hk + 1.0e-6, "final LP bound " + std::to_string(prev)
    + " does not exceed the tour optimum " + std::to_string(hk));

  /* Same for the MIP over those rows. The loop stops at connectivity, which in
   * general leaves some SECs unstated, so this is a relaxation and may come out
   * below hk with subtours -- what must not happen is coming out *above*. */
  Cbc_solve(m);
  check(Cbc_isProvenOptimal(m) != 0, "MIP over the loop's rows: proven optimal");
  double loopObj = Cbc_getObjValue(m);
  check(loopObj <= hk + 1.0e-6, "MIP over the loop's rows stays a relaxation ("
    + std::to_string(loopObj) + " <= " + std::to_string(hk) + ")");

  /* Now state the rest of the family lazily, on top of the loop's rows: the
   * combination must land exactly on the optimum. */
  std::vector<int> idx(static_cast<size_t>(n) * n);
  std::vector<double> cf(static_cast<size_t>(n) * n, 1.0);
  for (int mask = 1; mask < (1 << n); ++mask) {
    int sz = 0;
    for (int b = 0; b < n; ++b)
      if (mask & (1 << b))
        ++sz;
    if (sz < 2 || sz > n - 1)
      continue;
    int k = 0;
    for (int i = 0; i < n; ++i) {
      if (!(mask & (1 << i)))
        continue;
      for (int j = 0; j < n; ++j)
        if (i != j && (mask & (1 << j)))
          idx[k++] = t.col(i, j);
    }
    Cbc_addLazyConstraint(m, k, &idx[0], &cf[0], 'L', sz - 1);
  }
  Cbc_solve(m);
  check(Cbc_isProvenOptimal(m) != 0,
    "MIP with the loop's rows plus lazy SECs: proven optimal");
  checkTourSolution(m, t, hk, "MIP with the loop's rows plus lazy SECs");
  Cbc_deleteModel(m);
}

/* ── Test 9: staged rows and columns reaching the solve ────────────────── */

void testStagedModelReachesSolve(const Tsp &t)
{
  printf("\n=== Test 9: buffered rows/columns reach the solve ===\n");

  /* A model built purely through the staging API and handed straight to
   * Cbc_solve(). Cbc dispatches to the LP path when the solver reports no
   * integers, so an unflushed solver routes this to the LP branch and returns
   * the relaxation (0.5) while claiming proven optimality. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setObjSense(m, 1.0);
    Cbc_setLogLevel(m, 0);
    Cbc_addCol(m, "x", 0.0, 10.0, 1.0, 1, 0, NULL, NULL);
    int idx[1] = { 0 };
    double cf[1] = { 2.0 };
    Cbc_addRow(m, "c", 1, idx, cf, 'G', 1.0);
    Cbc_solve(m);
    checkEqInt(Cbc_getNumIntegers(m), 1, "fresh model: integer count seen by the solve");
    /* LP bound is 0.5; only a MIP solve returns 1. */
    checkClose(Cbc_getObjValue(m), 1.0, 1.0e-6,
      "fresh model solved as a MIP, not as its LP relaxation");
    Cbc_deleteModel(m);
  }

  /* A row added after an earlier solve must constrain the next one. Forbidding
   * every arc of the optimal tour, in both directions, makes the model strictly
   * worse -- unless the row never reached it. */
  {
    double hk = heldKarp(t);
    Cbc_Model *m = buildTsp(t, SEC_ROWS);
    Cbc_solve(m);
    checkClose(Cbc_getObjValue(m), hk, 1.0e-6, "first solve reaches the tour optimum");
    const double *x = Cbc_getColSolution(m);
    std::vector<int> idx;
    std::vector<double> cf;
    for (int i = 0; i < t.n; ++i)
      for (int j = 0; j < t.n; ++j)
        if (x[t.col(i, j)] > 0.5) {
          idx.push_back(t.col(i, j));
          cf.push_back(1.0);
          idx.push_back(t.col(j, i));
          cf.push_back(1.0);
        }
    int rowsBefore = Cbc_getNumRows(m);
    Cbc_addRow(m, "forbidTour", static_cast<int>(idx.size()), &idx[0], &cf[0], 'L',
      t.n - 1.0);
    Cbc_solve(m);
    checkEqInt(Cbc_getNumRows(m), rowsBefore + 1, "row added after a solve is in the model");
    check(Cbc_isProvenOptimal(m) != 0, "solve after adding the row: proven optimal");
    check(Cbc_getObjValue(m) > hk + 1.0e-6,
      "the added row binds: second optimum " + std::to_string(Cbc_getObjValue(m))
        + " is strictly worse than " + std::to_string(hk));
    Cbc_deleteModel(m);
  }
}

/* ── Test 10: the whole solution pool respects the lazy constraints ────── */

void testSolutionPool(const Tsp &t)
{
  printf("\n=== Test 10: every pooled solution respects the lazy constraints — n=%d ===\n",
    t.n);
  Cbc_Model *m = buildTsp(t, SEC_LAZY);
  Cbc_solve(m);
  int nSol = Cbc_numberSavedSolutions(m);
  check(nSol >= 1, "solution pool is not empty (" + std::to_string(nSol) + " solution(s))");

  int infeasible = 0, notTour = 0, objMismatch = 0;
  for (int s = 0; s < nSol; ++s) {
    const double *sol = Cbc_savedSolution(m, s);
    double maxViolRow = 0.0, maxViolCol = 0.0;
    int rowIdx = -1, colIdx = -1;
    if (!Cbc_checkFeasibility(m, sol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
      ++infeasible;
      printf("    pool solution %d/%d: rowViol=%.3e (row %d) colViol=%.3e (col %d)\n",
        s + 1, nSol, maxViolRow, rowIdx, maxViolCol, colIdx);
      continue;
    }
    double cost = 0.0;
    int nFrac = 0;
    if (!isSingleTour(t, sol, cost, nFrac)) {
      ++notTour;
      printf("    pool solution %d/%d: not one Hamiltonian cycle (frac=%d)\n",
        s + 1, nSol, nFrac);
      continue;
    }
    if (fabs(cost - Cbc_savedSolutionObj(m, s)) > 1.0e-6) {
      ++objMismatch;
      printf("    pool solution %d/%d: recomputed cost %.6g vs reported %.6g\n",
        s + 1, nSol, cost, Cbc_savedSolutionObj(m, s));
    }
  }
  checkEqInt(infeasible, 0, "pooled solutions failing Cbc_checkFeasibility");
  checkEqInt(notTour, 0, "pooled solutions that are not a single Hamiltonian cycle");
  checkEqInt(objMismatch, 0, "pooled solutions whose reported objective disagrees "
    "with the recomputed tour cost");
  Cbc_deleteModel(m);
}

/* ── Test 11: a solution offered from outside the search ───────────────── */

/* A cut generator that never cuts. Registered with atSolution=1 it still marks
 * the model as carrying an at-solution generator, which is all these cases
 * need -- the point is the screening path, not the cuts. */
void noopCutCallback(void * /*osiSolver*/, void * /*osiCuts*/, void * /*appData*/,
  int /*level*/, int /*pass*/)
{
}

/* Check the published point, not just the published number.
 *
 * Both failures these cases guard against are wrong answers, but they look
 * nothing alike. Screening the wrong vector rejects a good solution outright,
 * so the model claims proven infeasibility at 1.8e308 -- any objective check
 * sees it. Losing the integer fixing before the final confirming resolve
 * instead publishes the *root LP relaxation*: a plausible-looking objective, at
 * a point that is either fractional or, once INT_PARAM_ROUND_INT_VARS has
 * rounded it to look integral, violates the very row it was solved under. That
 * one is only caught by an objective check because the relaxation bound happens
 * to differ from the optimum on this model -- so verify the solution too, which
 * is what a caller would actually use. */
void checkPublishedPoint(Cbc_Model *m, double objCoef, double rowUb,
  const std::string &tag)
{
  const double *sol = Cbc_getColSolution(m);
  if (!sol) {
    check(false, tag + "solution published");
    return;
  }
  double activity = 0.0, cost = 0.0;
  bool integral = true;
  for (int j = 0; j < 2; ++j) {
    double v = sol[j];
    if (fabs(v - floor(v + 0.5)) > 1.0e-6)
      integral = false;
    activity += v;
    cost += objCoef * v;
  }
  check(integral, tag + "published solution is integral");
  check(activity <= rowUb + 1.0e-6, tag + "published solution satisfies the row");
  checkClose(cost, Cbc_getObjValue(m), 1.0e-6,
    tag + "objective recomputed from the published solution");
}

/* Four tiny models covering what happens to a solution that reaches the model
 * from outside the search -- an explicit MIPStart, or the previous answer
 * replayed on a re-solve -- when the model screens at solutions.
 *
 * Such a candidate is screened by re-running the at-solution generators. Those
 * read the solver's current column solution while the candidate lives in a
 * separate array, so the screen has to install the candidate first; otherwise
 * it judges whatever the solver happens to hold and rejects a perfectly good
 * solution the moment that leftover vector is fractional. Accepting it is then
 * the other half of the job: the search confirms the winner by fixing its
 * integer variables and re-optimising the continuous ones, and that fixing has
 * to survive to the resolve or the published answer is the relaxation.
 *
 * The TSP tests above cannot see either failure. There the leftover vector is
 * the previous solve's optimal tour -- integral, and satisfying every SEC -- so
 * reading the wrong array happens to give the right verdict, and the SEC-closed
 * relaxation is tight enough at the optimum to hide a relaxed resolve. What is
 * needed is a model whose leftover vector is fractional at screen time and
 * whose relaxation is strictly loose; two binaries under x + y <= 1.5 do both,
 * since the LP relaxation sits at x = y = 0.75.
 *
 * Kept in both formulations (an explicit lazy constraint, and an at-solution
 * cut callback with no lazy constraint at all) because the two arrive at the
 * screen by different routes, and in the maximize sense because the cutoff a
 * wrongly-accepted or wrongly-rejected solution installs is negated there. */
void testScreenReadsCandidate()
{
  printf("\n=== Test 11: a solution offered from outside the search — "
    "screened against the candidate, confirmed at its own bounds ===\n");

  /* (a) explicit lazy constraint; optimum is all-zero, so nothing is carried
   *     over as a MIPStart and each solve must stand on its own. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setLogLevel(m, 0);
    Cbc_setObjSense(m, -1.0); /* maximize */
    Cbc_addCol(m, "x", 0.0, 1.0, -1.0, 1, 0, NULL, NULL);
    Cbc_addCol(m, "y", 0.0, 1.0, -1.0, 1, 0, NULL, NULL);
    int idx[2] = { 0, 1 };
    double cf[2] = { 1.0, 1.0 };
    Cbc_addRow(m, "c", 2, idx, cf, 'L', 1.5);
    Cbc_addLazyConstraint(m, 2, idx, cf, 'L', 1.0);
    for (int s = 1; s <= 3; ++s) {
      Cbc_solve(m);
      std::string tag = "max + lazy, all-zero optimum, solve "
        + std::to_string(s) + ": ";
      check(Cbc_isProvenOptimal(m) != 0, tag + "proven optimal");
      check(Cbc_isProvenInfeasible(m) == 0, tag + "not reported infeasible");
      checkClose(Cbc_getObjValue(m), 0.0, 1.0e-6, tag + "objective");
      checkPublishedPoint(m, -1.0, 1.5, tag);
    }
    Cbc_deleteModel(m);
  }

  /* (b) at-solution cut callback and *no* lazy constraint. The model is still
   *     flagged as screening at solutions, so the same path runs -- and the
   *     leftover vector at screen time is the fractional LP relaxation. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setLogLevel(m, 0);
    Cbc_setObjSense(m, -1.0); /* maximize */
    Cbc_addCol(m, "x", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    Cbc_addCol(m, "y", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    int idx[2] = { 0, 1 };
    double cf[2] = { 1.0, 1.0 };
    Cbc_addRow(m, "c", 2, idx, cf, 'L', 1.5);
    Cbc_addCutCallback(m, noopCutCallback, "noop", NULL, 1, 1);
    for (int s = 1; s <= 3; ++s) {
      Cbc_solve(m);
      std::string tag = "max + atSolution cut callback, no lazy constraint, solve "
        + std::to_string(s) + ": ";
      check(Cbc_isProvenOptimal(m) != 0, tag + "proven optimal");
      check(Cbc_isProvenInfeasible(m) == 0, tag + "not reported infeasible");
      checkClose(Cbc_getObjValue(m), 1.0, 1.0e-6, tag + "objective");
      checkPublishedPoint(m, 1.0, 1.5, tag);
    }
    Cbc_deleteModel(m);
  }

  /* (c) the same callback registered with atSolution=0 -- an ordinary cut
   *     generator, which must not put the model on the screening path at all.
   *     This is the control: it passes either way, and its job is to show that
   *     (b) fails because of the at-solution flag and nothing else. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setLogLevel(m, 0);
    Cbc_setObjSense(m, -1.0); /* maximize */
    Cbc_addCol(m, "x", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    Cbc_addCol(m, "y", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    int idx[2] = { 0, 1 };
    double cf[2] = { 1.0, 1.0 };
    Cbc_addRow(m, "c", 2, idx, cf, 'L', 1.5);
    Cbc_addCutCallback(m, noopCutCallback, "noop", NULL, 1, 0);
    for (int s = 1; s <= 3; ++s) {
      Cbc_solve(m);
      std::string tag = "max + ordinary cut callback (control), solve "
        + std::to_string(s) + ": ";
      check(Cbc_isProvenOptimal(m) != 0, tag + "proven optimal");
      checkClose(Cbc_getObjValue(m), 1.0, 1.0e-6, tag + "objective");
      checkPublishedPoint(m, 1.0, 1.5, tag);
    }
    Cbc_deleteModel(m);
  }

  /* (d) the same model on its *first* solve, with the candidate supplied
   *     explicitly instead of arriving as the replayed previous answer. This is
   *     the minimal form of the whole situation: one Cbc_solve() on a fresh
   *     model, a partial MIPStart, and an at-solution generator. It removes any
   *     doubt that the earlier cases depend on state left over from a preceding
   *     solve. */
  {
    Cbc_Model *m = Cbc_newModel();
    Cbc_setLogLevel(m, 0);
    Cbc_setObjSense(m, -1.0); /* maximize */
    Cbc_addCol(m, "x", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    Cbc_addCol(m, "y", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);
    int idx[2] = { 0, 1 };
    double cf[2] = { 1.0, 1.0 };
    Cbc_addRow(m, "c", 2, idx, cf, 'L', 1.5);
    Cbc_addCutCallback(m, noopCutCallback, "noop", NULL, 1, 1);
    int msIdx[1] = { 1 };
    double msVal[1] = { 1.0 };
    Cbc_setMIPStartI(m, 1, msIdx, msVal); /* y = 1, x left to the search */
    Cbc_solve(m);
    std::string tag = "max + atSolution cut callback + explicit MIPStart: ";
    check(Cbc_isProvenOptimal(m) != 0, tag + "proven optimal");
    check(Cbc_isProvenInfeasible(m) == 0, tag + "not reported infeasible");
    checkClose(Cbc_getObjValue(m), 1.0, 1.0e-6, tag + "objective");
    checkPublishedPoint(m, 1.0, 1.5, tag);
    Cbc_deleteModel(m);
  }
}

} /* anonymous namespace */

int main(int argc, char **argv)
{
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "--help" || arg == "-h") {
      printf(
        "lazy-constraint-test — tests for Cbc's C-interface lazy constraints\n"
        "(Cbc_addLazyConstraint and Cbc_addCutCallback with atSolution=1),\n"
        "driven by the TSP with subtour elimination.\n"
        "\n"
        "Usage: lazy-constraint-test [--help]\n"
        "\n"
        "No data files are needed: instances come from a fixed LCG and every\n"
        "reference optimum is computed in-process three independent ways --\n"
        "Held-Karp dynamic programming, the same MIP with the SECs as ordinary\n"
        "rows, and the same MIP with them as lazy constraints -- which must\n"
        "agree, while the SEC-free assignment relaxation must come out strictly\n"
        "below. A lazy constraint that is silently ignored still yields an\n"
        "optimal-looking answer to a *relaxation*, which is what that gap\n"
        "detects. Checks cover:\n"
        "\n"
        "  1. lazy SECs binding exactly as row SECs do, with the published\n"
        "     solution walked as a successor array (an objective check alone\n"
        "     would pass on two subtours priced as a tour);\n"
        "  2. the LP relaxation staying untouched by them;\n"
        "  3. repeated Cbc_solve() on one model, where the previous answer comes\n"
        "     back as a MIPStart -- screened before the model is flagged as\n"
        "     carrying lazy constraints;\n"
        "  4. lazy constraints registered between solves;\n"
        "  5. Cbc_clone() carrying them, as a deep copy;\n"
        "  6. a MIPStart that violates one being rejected rather than trusted\n"
        "     (both objective senses), and a lazy-feasible one being kept;\n"
        "  7. the same family expressed through a cut callback at integer\n"
        "     solutions, separating on demand;\n"
        "  8. the iterative LP cut loop (solveLinearProgram/addRow/resolve);\n"
        "  9. staged rows and columns reaching the solve at all;\n"
        " 10. every solution in the pool respecting the lazy constraints, not\n"
        "     just the best one;\n"
        " 11. a solution offered from outside the search being screened against\n"
        "     the *candidate* rather than whatever vector the solver happens to\n"
        "     hold, and then confirmed at its own fixed bounds -- on tiny\n"
        "     maximize models whose leftover vector is fractional and whose\n"
        "     relaxation is strictly loose, where the first failure turns a\n"
        "     solved model into a proven-infeasible one and the second publishes\n"
        "     the root LP relaxation as if it were the optimum.\n"
        "\n"
        "Exit code: 0 = all checks passed, 1 = one or more checks failed.\n");
      return 0;
    }
    fprintf(stderr, "lazy-constraint-test: unrecognized argument '%s' "
      "(try --help)\n", arg.c_str());
    return 1;
  }

  /* n <= 12 keeps both the full SEC enumeration (2^n subsets) and the exact
   * Held-Karp DP cheap; seed 42 throughout, so every run is identical. */
  const Tsp t8(8, 42), t10(10, 42), t12(12, 42);

  testEnforcement(t8);
  testEnforcement(t10);
  testEnforcement(t12);

  testLpRelaxationUnaffected(t8);
  testLpRelaxationUnaffected(t12);

  testRepeatedSolve(t8);
  testRepeatedSolve(t10);

  testAddAfterSolve(t8);
  testAddAfterSolve(t10);

  testClone(t8);

  testMipStart(t8);
  testMipStart(t10);

  testCutCallbackAtSolution(t8);
  testCutCallbackAtSolution(t10);

  testLpCutLoop(t8);
  testLpCutLoop(t12);

  testStagedModelReachesSolve(t8);

  testSolutionPool(t8);

  testScreenReadsCandidate();

  printf("\n=== Summary ===\n");
  if (g_failures == 0) {
    printf("RESULT: OK (all checks passed)\n");
    return 0;
  }
  printf("RESULT: FAILED (%d check(s) failed)\n", g_failures);
  return 1;
}
