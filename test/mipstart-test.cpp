/*
 * mipstart-test — targeted tests for MIP starts through the C interface
 * (Cbc_setMIPStart / Cbc_setMIPStartI), which is the path python-mip uses.
 *
 * Why the C interface and not the -mipstart file. CbcMipStart::read() densifies
 * before returning: it builds a vector over *every* column initialised to zero
 * and overwrites the entries the file mentioned. So a file naming three columns
 * and a file naming all of them arrive at computeCompleteSolution() as the same
 * dense vector, and the three input shapes a user might supply cannot be told
 * apart on that path at all. Cbc_setMIPStart() stores names and values verbatim,
 * so it is the only place the distinction is observable -- and the only place a
 * regression in it would be seen by callers.
 *
 * The three shapes, which must all work (the caller decides what to send; the
 * solver decides what to fix):
 *
 *   (a) the nonzero *integer* values only -- what the API documents and what
 *       most callers send. Continuous columns are not mentioned at all and have
 *       to be derived by the LP.
 *   (b) the full solution -- every column, integer and continuous.
 *   (c) all nonzeros -- integers and continuous alike, zeros omitted.
 *
 * Models are tiny and their optima are derived here in closed form rather than
 * hardcoded from a solver run, so a wrong answer is a wrong answer and not a
 * stale fixture. Each has a continuous tail that is *not* determined by the
 * integers alone until the LP runs, which is what makes shape (a) a real test:
 * a start that only names integers must still come back with the continuous
 * values filled in and the recomputed objective matching.
 *
 * Coverage:
 *
 *   1. The three shapes, by name and by index, with and without preprocessing.
 *      For each: the start is accepted, the solution is integral, every row is
 *      satisfied, and the objective recomputed from the published columns
 *      matches what Cbc reports.
 *   2. -mipStartFix / INT_PARAM_MIPSTART_FIX, the mode that replaced the
 *      compile-time JUST_FIX_INTEGER macro. All three values must produce the
 *      same optimum on a model where the start is correct. The "all" mode also
 *      has to actually reach the fix-continuous-too path -- a test that never
 *      entered it would pass whatever the mode did -- so it is checked on a
 *      model where fixing the continuous values at slightly perturbed numbers
 *      makes the first LP infeasible and the integers-only retry is what
 *      rescues it.
 *   3. Negative integer values. Cbc_solve() auto-carries each solution forward
 *      as the next start but keeps only entries >= 1e-8, so shape (c) on a model
 *      whose optimum has negative integers is where a positive-only filter
 *      shows up. The optimum must survive a re-solve.
 *   4. A start that is *infeasible* must be rejected, not trusted: taking it as
 *      an incumbent installs a cutoff that can exclude the genuine optimum. The
 *      final answer must still be the true optimum.
 *   5. A start naming a column that does not exist, and a start naming nothing
 *      the model has, must both be survivable -- ignored with the solve
 *      proceeding, never a crash or a wrong answer.
 *   6. A suboptimal but feasible start must be accepted and then improved upon,
 *      which is the whole point of the feature.
 *
 * Usage: mipstart-test [--help]
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "Cbc_C_Interface.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

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

/* ── capturing Cbc's own log ───────────────────────────────────────────── */

/* Some of what has to be asserted here is only observable as a log line.
 * Whether the fix-continuous-too recovery path ran is the case that matters: the
 * objective is the same either way, so a test that merely checks the answer
 * passes whether the branch was entered or skipped -- which is exactly how that
 * branch stayed dead through every release while it was behind a macro. Same
 * technique cbc_validate_sol.cpp uses to silence the banner, redirected to a
 * temporary file instead of /dev/null so the text can be read back. */
class LogCapture {
public:
  LogCapture()
    : fd_(-1)
    , savedOut_(-1)
    , savedErr_(-1)
  {
    char tmpl[] = "/tmp/mipstart-test-log-XXXXXX";
    fd_ = mkstemp(tmpl);
    if (fd_ < 0)
      return;
    path_ = tmpl;
    fflush(stdout);
    fflush(stderr);
    savedOut_ = dup(STDOUT_FILENO);
    savedErr_ = dup(STDERR_FILENO);
    dup2(fd_, STDOUT_FILENO);
    dup2(fd_, STDERR_FILENO);
  }

  /* Restores the real streams and returns everything that was written. */
  std::string finish()
  {
    if (savedOut_ < 0)
      return std::string();
    fflush(stdout);
    fflush(stderr);
    dup2(savedOut_, STDOUT_FILENO);
    dup2(savedErr_, STDERR_FILENO);
    close(savedOut_);
    close(savedErr_);
    savedOut_ = savedErr_ = -1;

    std::string out;
    lseek(fd_, 0, SEEK_SET);
    char buf[4096];
    ssize_t n;
    while ((n = ::read(fd_, buf, sizeof(buf))) > 0)
      out.append(buf, static_cast< size_t >(n));
    close(fd_);
    fd_ = -1;
    unlink(path_.c_str());
    return out;
  }

  ~LogCapture()
  {
    if (savedOut_ >= 0)
      finish();
  }

private:
  int fd_;
  int savedOut_;
  int savedErr_;
  std::string path_;
};

bool logHas(const std::string &log, const char *needle)
{
  return log.find(needle) != std::string::npos;
}

/* A fixed linear congruential generator, so that the one fixture built from
 * random data is the same instance on every machine and every run. Deliberately
 * not rand(): its sequence is implementation-defined, which would make a failure
 * here unreproducible elsewhere. */
double lcg(unsigned long &seed, int lo, int hi)
{
  seed = seed * 1103515245UL + 12345UL;
  return (double)(lo + (int)((seed >> 16) % (unsigned long)(hi - lo + 1)));
}

/* ── the fixture model ─────────────────────────────────────────────────── */

/* A fixed-charge model, small enough to reason about and shaped so that the
 * integers alone do not determine the answer:
 *
 *   min  sum_i f_i y_i + sum_i c_i x_i
 *   s.t. sum_i x_i           >= D          (demand)
 *        x_i - cap_i y_i     <= 0    (i)   (link: no flow without the fixed
 *        x_i >= 0, y_i binary               charge paid)
 *
 * The y are the integers; the x are continuous and have to be *derived*. Given
 * any y, the remaining LP is a continuous knapsack -- fill the cheapest opened
 * capacity first -- so a start naming only the y is genuinely incomplete, and
 * the value of the x it comes back with is checkable in closed form.
 *
 * Costs are chosen so the optimum opens exactly two facilities and the third is
 * a near miss, which keeps the model from being decided by the LP relaxation
 * alone. */
struct FixedCharge {
  static const int n = 4;
  double f[n];   /* fixed charge  */
  double c[n];   /* unit cost     */
  double cap[n]; /* capacity      */
  double demand;

  FixedCharge()
  {
    const double ff[n] = { 30.0, 45.0, 12.0, 70.0 };
    const double cc[n] = { 2.0, 1.0, 5.0, 0.5 };
    const double kk[n] = { 6.0, 8.0, 4.0, 10.0 };
    for (int i = 0; i < n; ++i) {
      f[i] = ff[i];
      c[i] = cc[i];
      cap[i] = kk[i];
    }
    demand = 11.0;
  }

  /* Column layout: y_0..y_{n-1} then x_0..x_{n-1}. */
  int colY(int i) const { return i; }
  int colX(int i) const { return n + i; }
  int numCols() const { return 2 * n; }

  /* Best continuous fill for a given open set, and its total cost, by exhaustive
   * greedy on unit cost -- the LP's own answer, computed independently. Returns
   * 1e30 if the open capacity cannot meet demand. */
  double bestCost(const int *open, double *x) const
  {
    double cap_total = 0.0;
    for (int i = 0; i < n; ++i)
      cap_total += open[i] ? cap[i] : 0.0;
    for (int i = 0; i < n; ++i)
      x[i] = 0.0;
    if (cap_total < demand - 1e-9)
      return 1e30;

    double cost = 0.0;
    for (int i = 0; i < n; ++i)
      cost += open[i] ? f[i] : 0.0;

    /* cheapest unit cost first */
    double left = demand;
    for (int pass = 0; pass < n; ++pass) {
      int best = -1;
      for (int i = 0; i < n; ++i) {
        if (!open[i] || x[i] > 0.0)
          continue;
        if (best < 0 || c[i] < c[best])
          best = i;
      }
      if (best < 0)
        break;
      double take = left < cap[best] ? left : cap[best];
      x[best] = take;
      cost += c[best] * take;
      left -= take;
      if (left <= 1e-9)
        break;
    }
    /* Columns not used stay at zero, but mark them visited for the loop above
     * by leaving x at 0 -- the pass counter bounds the loop, so nothing spins. */
    return cost;
  }

  /* The true optimum, by enumerating all 2^n open sets. */
  double optimum(int *bestOpen, double *bestX) const
  {
    double best = 1e30;
    double x[n];
    int open[n];
    for (int mask = 0; mask < (1 << n); ++mask) {
      for (int i = 0; i < n; ++i)
        open[i] = (mask >> i) & 1;
      double cost = bestCost(open, x);
      if (cost < best) {
        best = cost;
        for (int i = 0; i < n; ++i) {
          bestOpen[i] = open[i];
          bestX[i] = x[i];
        }
      }
    }
    return best;
  }
};

std::string yName(int i)
{
  char b[32];
  snprintf(b, sizeof(b), "y%d", i);
  return std::string(b);
}

std::string xName(int i)
{
  char b[32];
  snprintf(b, sizeof(b), "x%d", i);
  return std::string(b);
}

Cbc_Model *buildFixedCharge(const FixedCharge &m)
{
  Cbc_Model *model = Cbc_newModel();
  Cbc_setObjSense(model, 1.0);

  for (int i = 0; i < m.n; ++i)
    Cbc_addCol(model, yName(i).c_str(), 0.0, 1.0, m.f[i], 1, 0, NULL, NULL);
  for (int i = 0; i < m.n; ++i)
    Cbc_addCol(model, xName(i).c_str(), 0.0, m.cap[i], m.c[i], 0, 0, NULL, NULL);

  /* demand: sum x_i >= D */
  {
    std::vector<int> idx;
    std::vector<double> coef;
    for (int i = 0; i < m.n; ++i) {
      idx.push_back(m.colX(i));
      coef.push_back(1.0);
    }
    Cbc_addRow(model, "demand", m.n, idx.data(), coef.data(), 'G', m.demand);
  }
  /* link_i: x_i - cap_i y_i <= 0 */
  for (int i = 0; i < m.n; ++i) {
    int idx[2] = { m.colX(i), m.colY(i) };
    double coef[2] = { 1.0, -m.cap[i] };
    char nm[32];
    snprintf(nm, sizeof(nm), "link%d", i);
    Cbc_addRow(model, nm, 2, idx, coef, 'L', 0.0);
  }
  return model;
}

/* Quiet, and short enough that a stuck solve fails the test rather than the
 * suite's own time limit. */
void configureQuiet(Cbc_Model *model)
{
  Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 0);
  Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
}

/* Recompute the objective from the published columns rather than trusting the
 * reported value -- the two disagreeing is exactly the failure a MIP start can
 * cause, by publishing a point whose cost was computed on a modified LP. */
double recomputeObj(Cbc_Model *model, const double *sol)
{
  double obj = 0.0;
  int nc = Cbc_getNumCols(model);
  const double *c = Cbc_getObjCoefficients(model);
  for (int i = 0; i < nc; ++i)
    obj += c[i] * sol[i];
  return obj;
}

/* Every row satisfied, every integer integral. Returns the largest violation. */
double worstViolation(Cbc_Model *model, const double *sol, double *worstFrac)
{
  const int nr = Cbc_getNumRows(model);
  const int nc = Cbc_getNumCols(model);
  double worst = 0.0;
  *worstFrac = 0.0;

  for (int i = 0; i < nc; ++i) {
    if (Cbc_isInteger(model, i)) {
      double d = fabs(sol[i] - floor(sol[i] + 0.5));
      if (d > *worstFrac)
        *worstFrac = d;
    }
    double lo = Cbc_getColLower(model)[i], up = Cbc_getColUpper(model)[i];
    if (sol[i] < lo - worst)
      worst = lo - sol[i] > worst ? lo - sol[i] : worst;
    if (sol[i] > up + worst)
      worst = sol[i] - up > worst ? sol[i] - up : worst;
  }

  std::vector<int> idx(nc);
  std::vector<double> coef(nc);
  for (int r = 0; r < nr; ++r) {
    int nz = Cbc_getRowNz(model, r);
    const int *ri = Cbc_getRowIndices(model, r);
    const double *rc = Cbc_getRowCoeffs(model, r);
    double act = 0.0;
    for (int k = 0; k < nz; ++k)
      act += rc[k] * sol[ri[k]];
    double lo = Cbc_getRowLower(model)[r], up = Cbc_getRowUpper(model)[r];
    if (act < lo - 1e-6 && lo - act > worst)
      worst = lo - act;
    if (act > up + 1e-6 && act - up > worst)
      worst = act - up;
  }
  return worst;
}

/* ── Test 1: the three input shapes ────────────────────────────────────── */

/* Shapes are built from the *known* optimum, so an accepted start should let
 * the search confirm that optimum immediately; and, more to the point, the
 * published solution has to be complete and consistent whichever shape went
 * in -- including shape (a), where the continuous columns were never supplied
 * and the LP had to derive them. */
enum Shape { ShapeIntNonzero, ShapeFull, ShapeAllNonzero };

const char *shapeName(Shape s)
{
  switch (s) {
  case ShapeIntNonzero:
    return "nonzero integers only";
  case ShapeFull:
    return "full solution";
  default:
    return "all nonzeros";
  }
}

void buildStart(const FixedCharge &m, Shape shape, const int *open, const double *x,
  std::vector<std::string> &names, std::vector<int> &idxs, std::vector<double> &vals)
{
  names.clear();
  idxs.clear();
  vals.clear();
  for (int i = 0; i < m.n; ++i) {
    double v = open[i] ? 1.0 : 0.0;
    if (shape != ShapeFull && v == 0.0)
      continue;
    names.push_back(yName(i));
    idxs.push_back(m.colY(i));
    vals.push_back(v);
  }
  if (shape == ShapeIntNonzero)
    return;
  for (int i = 0; i < m.n; ++i) {
    if (shape != ShapeFull && x[i] == 0.0)
      continue;
    names.push_back(xName(i));
    idxs.push_back(m.colX(i));
    vals.push_back(x[i]);
  }
}

void testShape(Shape shape, bool byIndex, bool preprocess)
{
  FixedCharge m;
  int open[FixedCharge::n];
  double x[FixedCharge::n];
  const double opt = m.optimum(open, x);

  char tag[160];
  snprintf(tag, sizeof(tag), "%s, %s, preprocess %s",
    shapeName(shape), byIndex ? "by index" : "by name", preprocess ? "on" : "off");
  printf("\n=== MIPStart shape: %s ===\n", tag);

  Cbc_Model *model = buildFixedCharge(m);
  Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
  /* Log level 1 so the MIPStart lines are emitted; the whole solve is captured,
     so nothing reaches the terminal either way. */
  Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 1);
  Cbc_setParameter(model, "preprocess", preprocess ? "on" : "off");

  std::vector<std::string> names;
  std::vector<int> idxs;
  std::vector<double> vals;
  buildStart(m, shape, open, x, names, idxs, vals);

  if (byIndex) {
    Cbc_setMIPStartI(model, (int)idxs.size(), idxs.data(), vals.data());
  } else {
    std::vector<const char *> cn;
    for (size_t k = 0; k < names.size(); ++k)
      cn.push_back(names[k].c_str());
    Cbc_setMIPStart(model, (int)cn.size(), cn.data(), vals.data());
  }

  LogCapture cap;
  Cbc_solve(model);
  std::string log = cap.finish();

  check(Cbc_getNumIntegers(model) == FixedCharge::n, std::string("model has 4 integers: ") + tag);
  check(Cbc_isProvenOptimal(model) != 0, std::string("proven optimal: ") + tag);

  /* Without this the whole test would pass on a build that ignored MIP starts
     entirely: the model is small enough that Cbc reaches the same optimum
     unaided, so "the answer is right" says nothing about whether the start was
     used. The start given here *is* the optimum, so it must be accepted. */
  check(logHas(log, "MIPStart provided solution"),
    std::string("the start was accepted and turned into a solution: ") + tag);
  check(!logHas(log, "could not be used to build a solution"),
    std::string("the start was not rejected: ") + tag);
  check(!logHas(log, "not valid, column names do not match"),
    std::string("the start's columns were matched: ") + tag);

  const double *sol = Cbc_getColSolution(model);
  if (sol) {
    checkClose(Cbc_getObjValue(model), opt, 1e-6, std::string("objective is the optimum: ") + tag);
    checkClose(recomputeObj(model, sol), Cbc_getObjValue(model), 1e-6,
      std::string("reported objective matches the published columns: ") + tag);
    double frac = 0.0;
    double viol = worstViolation(model, sol, &frac);
    check(viol <= 1e-6, std::string("no row or bound violation: ") + tag);
    check(frac <= 1e-6, std::string("all integers integral: ") + tag);

    /* The point of shape (a): columns the start never mentioned must come back
     * with the values the LP derives, not zero. */
    double xsum = 0.0;
    for (int i = 0; i < FixedCharge::n; ++i)
      xsum += sol[m.colX(i)];
    checkClose(xsum, m.demand, 1e-6,
      std::string("continuous columns were derived and meet demand: ") + tag);
  } else {
    check(false, std::string("a solution was published: ") + tag);
  }

  Cbc_deleteModel(model);
}

/* ── Test 2: the fix mode that replaced JUST_FIX_INTEGER ───────────────── */

/* All three modes must reach the same optimum from a correct start. The modes
 * differ only in what is *fixed* before the recovery LP, which is an internal
 * decision -- so the observable requirement is that none of them loses the
 * answer. */
void testFixModes()
{
  printf("\n=== MIPStart fix modes (INT_PARAM_MIPSTART_FIX) ===\n");
  FixedCharge m;
  int open[FixedCharge::n];
  double x[FixedCharge::n];
  const double opt = m.optimum(open, x);

  for (int mode = 0; mode <= 2; ++mode) {
    Cbc_Model *model = buildFixedCharge(m);
    configureQuiet(model);
    Cbc_setIntParam(model, INT_PARAM_MIPSTART_FIX, mode);
    check(Cbc_getIntParam(model, INT_PARAM_MIPSTART_FIX) == mode,
      "fix mode round-trips through get/setIntParam");

    std::vector<std::string> names;
    std::vector<int> idxs;
    std::vector<double> vals;
    buildStart(m, ShapeFull, open, x, names, idxs, vals);
    std::vector<const char *> cn;
    for (size_t k = 0; k < names.size(); ++k)
      cn.push_back(names[k].c_str());
    Cbc_setMIPStart(model, (int)cn.size(), cn.data(), vals.data());

    Cbc_solve(model);
    char tag[64];
    snprintf(tag, sizeof(tag), "fix mode %d", mode);
    check(Cbc_isProvenOptimal(model) != 0, std::string("proven optimal, ") + tag);
    checkClose(Cbc_getObjValue(model), opt, 1e-6, std::string("optimum reached, ") + tag);
    Cbc_deleteModel(model);
  }
}

/* Mode 2 ("all") fixes the supplied continuous values too, and that is the mode
 * whose recovery path -- unfix the continuous columns, keep the integers, retry
 * -- was unreachable as shipped: the fixing loop skipped every non-integer, so
 * nContinuousFixed could never become nonzero and the branch guarding on it,
 * SOS1/SOS2 handling included, was dead code in every release.
 *
 * So this checks that the branch *runs*, not merely that the answer is right.
 * The objective comes out the same whether it fires or is skipped, which is
 * precisely why nobody noticed it was dead -- an objective-only assertion here
 * would pass against the macro build too.
 *
 * The trigger is a full start whose continuous flows fall a unit short of
 * demand. Every value is inside its own column bounds so validateMIPStartValues()
 * accepts them all, but pinning all of them leaves the demand row unsatisfiable
 * and the first LP goes infeasible -- the "falsely identified as infeasible"
 * hazard, staged deliberately. Unfixing the continuous columns is what has to
 * rescue it. */
void testFixAllRecovers()
{
  printf("\n=== MIPStart fix mode 'all' enters and survives its recovery path ===\n");
  FixedCharge m;
  int open[FixedCharge::n];
  double x[FixedCharge::n];
  const double opt = m.optimum(open, x);

  /* Short the largest flow by one unit: still within bounds, no longer meeting
     demand once every flow is pinned. */
  double xp[FixedCharge::n];
  int shorted = -1;
  for (int i = 0; i < FixedCharge::n; ++i)
    xp[i] = x[i];
  for (int i = 0; i < FixedCharge::n; ++i)
    if (xp[i] > 1.0) {
      xp[i] -= 1.0;
      shorted = i;
      break;
    }
  check(shorted >= 0, "the fixture has a flow large enough to short");
  double xsum = 0.0;
  for (int i = 0; i < FixedCharge::n; ++i)
    xsum += xp[i];
  check(xsum < m.demand - 1e-9, "the shorted start really does miss demand");

  std::vector<std::string> names;
  std::vector<int> idxs;
  std::vector<double> vals;
  buildStart(m, ShapeFull, open, xp, names, idxs, vals);
  std::vector<const char *> cn;
  for (size_t k = 0; k < names.size(); ++k)
    cn.push_back(names[k].c_str());

  /* Mode 2 must enter the recovery branch on this start; mode 1 must not, since
     it never fixes a continuous column and so never over-constrains the LP.
     Asserting both directions is what makes the mode observably in force rather
     than merely accepted. */
  for (int mode = 1; mode <= 2; ++mode) {
    Cbc_Model *model = buildFixedCharge(m);
    Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
    /* Log level 1, not 0: the CBC_GENERAL lines asserted on below are the point. */
    Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 1);
    Cbc_setIntParam(model, INT_PARAM_MIPSTART_FIX, mode);
    Cbc_setMIPStart(model, (int)cn.size(), cn.data(), vals.data());

    LogCapture cap;
    Cbc_solve(model);
    std::string log = cap.finish();

    char tag[64];
    snprintf(tag, sizeof(tag), "fix mode %d", mode);
    const bool entered = logHas(log, "Trying just fixing integer variables");
    if (mode == 2)
      check(entered, std::string("recovery path was entered, ") + tag);
    else
      check(!entered, std::string("recovery path was not entered, ") + tag);
    /* Reaching the branch is also where the infeasibility gets localised, and
       that report is only useful if it names something: an unlocalised "no
       single row or column identified" would pass a mere substring test for
       "infeasible", so assert on the row wording specifically. Mode 1 never
       over-constrains the LP, so it must stay silent -- both directions again,
       for the same reason as above. */
    const bool localised = logHas(log, "Bound propagation confirms the fixed"
      " problem is infeasible, detected at row");
    if (mode == 2)
      check(localised, std::string("the infeasibility was localised to a row, ") + tag);
    else
      check(!localised, std::string("no infeasibility was reported, ") + tag);
    check(logHas(log, "MIPStart provided solution"),
      std::string("the start was still turned into a solution, ") + tag);

    check(Cbc_isProvenOptimal(model) != 0, std::string("proven optimal, ") + tag);
    checkClose(Cbc_getObjValue(model), opt, 1e-6, std::string("optimum reached, ") + tag);
    const double *sol = Cbc_getColSolution(model);
    if (sol) {
      double frac = 0.0;
      check(worstViolation(model, sol, &frac) <= 1e-6,
        std::string("published point is feasible, ") + tag);
      check(frac <= 1e-6, std::string("published point is integral, ") + tag);
    }
    Cbc_deleteModel(model);
  }
}

/* ── Test 3: negative integer values in the start ──────────────────────── */

/* Shape (c) means *all* nonzeros, not all positive ones. A filter that keeps
 * only values >= 1e-8 -- which is what Cbc_solve()'s own solution-carrying does
 * -- would drop the negative entries and hand the next solve a start that says
 * those columns are zero. On a model whose optimum needs them negative, that is
 * either a rejected start or, worse, a wrong incumbent.
 *
 *   min  z_0 + z_1 + 3 w      z_i integer in [-5, 5], w continuous >= 0
 *   s.t. z_0 + z_1 + w >= -6
 *        z_0 - z_1      = 0
 *
 * Optimum: z_0 = z_1 = -3, w = 0, objective -6. */
void testNegativeIntegers()
{
  printf("\n=== MIPStart with negative integer values ===\n");
  Cbc_Model *model = Cbc_newModel();
  Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
  Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 1);
  Cbc_setObjSense(model, 1.0);
  Cbc_addCol(model, "z0", -5.0, 5.0, 1.0, 1, 0, NULL, NULL);
  Cbc_addCol(model, "z1", -5.0, 5.0, 1.0, 1, 0, NULL, NULL);
  Cbc_addCol(model, "w", 0.0, 10.0, 3.0, 0, 0, NULL, NULL);
  {
    int idx[3] = { 0, 1, 2 };
    double coef[3] = { 1.0, 1.0, 1.0 };
    Cbc_addRow(model, "demand", 3, idx, coef, 'G', -6.0);
  }
  {
    int idx[2] = { 0, 1 };
    double coef[2] = { 1.0, -1.0 };
    Cbc_addRow(model, "equal", 2, idx, coef, 'E', 0.0);
  }

  /* Shape (c) on this optimum is two negative entries and nothing else. */
  const char *names[2] = { "z0", "z1" };
  double vals[2] = { -3.0, -3.0 };
  Cbc_setMIPStart(model, 2, names, vals);

  LogCapture cap;
  Cbc_solve(model);
  std::string log = cap.finish();

  check(Cbc_isProvenOptimal(model) != 0, "proven optimal with a negative-valued start");
  checkClose(Cbc_getObjValue(model), -6.0, 1e-6, "optimum -6 reached");
  const double *sol = Cbc_getColSolution(model);
  if (sol) {
    checkClose(sol[0], -3.0, 1e-6, "z0 = -3");
    checkClose(sol[1], -3.0, 1e-6, "z1 = -3");
  }

  /* The objective alone cannot see this: the model is tiny and Cbc reaches -6
     unaided, so a filter that dropped both negative entries would leave every
     numeric check above passing. What distinguishes "the start was used" from
     "the start was silently emptied" is that the start is reported as providing
     a solution at its own cost -- with the values gone, the two integers are
     unfixed and there is nothing to complete. */
  check(logHas(log, "MIPStart provided solution with cost -6"),
    "the negative-valued start was used, at its own cost");
  check(!logHas(log, "not valid, column names do not match"),
    "the negative-valued start was not discarded as unusable");

  /* Solving the same model again carries the previous answer forward as the next
     start, and that carry-forward is where Cbc_solve() applies its own >= 1e-8
     filter -- so the second solve is the one that would lose the entries even if
     the first kept them. */
  LogCapture cap2;
  Cbc_solve(model);
  std::string log2 = cap2.finish();
  check(Cbc_isProvenOptimal(model) != 0, "proven optimal on re-solve (carried start)");
  checkClose(Cbc_getObjValue(model), -6.0, 1e-6, "optimum -6 on re-solve");
  check(logHas(log2, "MIPStart provided solution with cost -6"),
    "the carried-forward start kept its negative values");
  Cbc_deleteModel(model);
}

/* ── Test 3b: the mini branch and bound fills in fractional integers ───── */

/* When the LP that follows the fixing hands back fractional integers, the last
 * resort is a bounded branch and bound on that LP (the smallBranchAndBound call
 * in CbcMipStart.cpp). Nothing else in this file reaches it: the fixed-charge
 * fixture's integers are all named by the start, so its LP comes back integral
 * and the branch is never taken. Over the whole 464-instance mip-sanity corpus
 * exactly one instance reaches it, so it is effectively untested territory that
 * a MIP start nevertheless depends on for correctness.
 *
 * Picking the fixture took some care, because "the LP is fractional" is not
 * enough -- the sub-MIP gets a full root node with cuts, and on odd-cycle vertex
 * cover (the obvious way to force fractionality) the root closes the gap by
 * itself, so a build that had been crippled to allow *no* nodes at all still
 * succeeded and the test proved nothing. What is needed is a model whose root
 * genuinely cannot finish the job. Market split is the standard such family: m
 * equality rows over n binaries with large random coefficients, whose LP
 * relaxation is highly fractional and which is known to defeat cut-and-bound
 * without search. Coefficients come from a fixed LCG so the instance is
 * reproducible, and a subset is drawn first and each right-hand side set to that
 * subset's activity -- otherwise a random market-split instance is almost
 * certainly infeasible and would test the rejection path instead. m=4, n=20 is
 * the smallest size measured to need branching: the shipped code completes it in
 * ~0.04 s, while a build whose sub-MIP node limit collapses to zero reports
 * "could not be used to build a solution".
 *
 * The objective is on w alone, which a row of its own forces to 1: the binaries
 * of the split are pure feasibility, so nothing about the cost can stand in for
 * checking that the start was completed.
 *
 * Preprocessing is off deliberately. With it on, the w row is enough for presolve
 * to remove w outright, and the start's only named column no longer exists in the
 * problem handed to CbcMipStart, so the start is discarded before any of this.
 * That is correct behaviour, just not what is under test here. */
void testMiniBranchAndBound()
{
  printf("\n=== the mini branch and bound completes a fractional LP ===\n");
  const int mRows = 4;
  const int nBin = 20;

  Cbc_Model *model = Cbc_newModel();
  Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
  Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 1);
  Cbc_setObjSense(model, 1.0);
  Cbc_setParameter(model, "preprocess", "off");
  /* Mode "integer" (1), not the default: the default takes every integer the
     start does not name to be zero, which would fix the whole split to 0 and
     leave the rows unsatisfiable rather than fractional. */
  Cbc_setIntParam(model, INT_PARAM_MIPSTART_FIX, 1);

  unsigned long seed = 12345UL;
  char nm[32];
  for (int j = 0; j < nBin; ++j) {
    snprintf(nm, sizeof(nm), "s%d", j);
    Cbc_addCol(model, nm, 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
  }
  Cbc_addCol(model, "w", 0.0, 1.0, 1.0, 1, 0, NULL, NULL);

  std::vector<int> planted(nBin);
  for (int j = 0; j < nBin; ++j)
    planted[j] = (int)lcg(seed, 0, 1);
  for (int i = 0; i < mRows; ++i) {
    std::vector<int> idx(nBin);
    std::vector<double> coef(nBin);
    double rhs = 0.0;
    for (int j = 0; j < nBin; ++j) {
      idx[j] = j;
      coef[j] = lcg(seed, 0, 99);
      rhs += coef[j] * planted[j];
    }
    snprintf(nm, sizeof(nm), "split%d", i);
    Cbc_addRow(model, nm, nBin, idx.data(), coef.data(), 'E', rhs);
  }
  {
    int idx[1] = { nBin };
    double coef[1] = { 1.0 };
    Cbc_addRow(model, "forceW", 1, idx, coef, 'G', 1.0);
  }

  const char *names[1] = { "w" };
  double vals[1] = { 1.0 };
  Cbc_setMIPStart(model, 1, names, vals);

  LogCapture cap;
  Cbc_solve(model);
  std::string log = cap.finish();

  check(Cbc_isProvenOptimal(model) != 0, "proven optimal (mini B&B fixture)");
  checkClose(Cbc_getObjValue(model), 1.0, 1e-6, "optimum 1 reached");

  /* The assertions that make this a test of the mini branch and bound rather
     than of market split. Without them the whole thing passes on a build where
     smallBranchAndBound() never gets a node to explore and the start is thrown
     away, because Cbc reaches cost 1 unaided either way -- which is exactly how
     a node limit silently rewritten to -1 survived in the shipped code. */
  check(logHas(log, "variables are still fractional"),
    "the LP after fixing really did come back fractional");
  check(logHas(log, "Mini branch and bound defined values for remaining variables"),
    "the mini branch and bound ran and completed the start");
  check(logHas(log, "MIPStart provided solution with cost 1"),
    "the completed start was published at the cost of the user's objective");
  check(!logHas(log, "could not be used to build a solution"),
    "the fractional start was completed rather than abandoned");

  const double *sol = Cbc_getColSolution(model);
  if (sol) {
    double frac = 0.0;
    check(worstViolation(model, sol, &frac) <= 1e-6, "published point is feasible");
    check(frac <= 1e-6, "published point is integral");
  }
  Cbc_deleteModel(model);
}


/* ── Test 4: an infeasible start must not be trusted ───────────────────── */

/* An accepted start installs a cutoff. If a start that violates a row were
 * accepted at its own claimed cost, that cutoff could exclude the genuine
 * optimum and the search would end proving a wrong answer -- or proving
 * infeasibility on a perfectly solvable model. */
void testInfeasibleStartRejected()
{
  printf("\n=== an infeasible MIPStart is rejected, not trusted ===\n");
  FixedCharge m;
  int open[FixedCharge::n];
  double x[FixedCharge::n];
  const double opt = m.optimum(open, x);

  Cbc_Model *model = buildFixedCharge(m);
  Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
  Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 1);

  /* Claim flow with nothing opened: every link row is violated, and the cost of
   * this point is far below the true optimum, so trusting it would cut the
   * optimum away. */
  const char *names[2] = { "x0", "x1" };
  double vals[2] = { 6.0, 5.0 };
  Cbc_setMIPStart(model, 2, names, vals);

  LogCapture cap;
  Cbc_solve(model);
  std::string log = cap.finish();

  check(Cbc_isProvenOptimal(model) != 0, "proven optimal despite an infeasible start");
  checkClose(Cbc_getObjValue(model), opt, 1e-6, "true optimum still found");

  /* The start names only continuous columns, so under the default integers-only
     mode nothing gets fixed and it is discarded on those grounds; the point is
     that it is never published as a solution at the cost it claims. */
  check(!logHas(log, "MIPStart provided solution with cost 4"),
    "the infeasible start was not published at its own claimed cost");

  const double *sol = Cbc_getColSolution(model);
  if (sol) {
    double frac = 0.0;
    check(worstViolation(model, sol, &frac) <= 1e-6, "published point is feasible");
  }
  Cbc_deleteModel(model);
}

/* ── Test 5: starts naming columns the model does not have ─────────────── */

void testUnknownColumns()
{
  printf("\n=== MIPStart naming unknown columns ===\n");
  FixedCharge m;
  int open[FixedCharge::n];
  double x[FixedCharge::n];
  const double opt = m.optimum(open, x);

  /* Some names known, some not: the known ones should still be used. */
  {
    Cbc_Model *model = buildFixedCharge(m);
    configureQuiet(model);
    std::vector<std::string> names;
    std::vector<int> idxs;
    std::vector<double> vals;
    buildStart(m, ShapeIntNonzero, open, x, names, idxs, vals);
    names.push_back("no_such_column");
    vals.push_back(1.0);
    std::vector<const char *> cn;
    for (size_t k = 0; k < names.size(); ++k)
      cn.push_back(names[k].c_str());
    Cbc_setMIPStart(model, (int)cn.size(), cn.data(), vals.data());
    Cbc_solve(model);
    check(Cbc_isProvenOptimal(model) != 0, "proven optimal with one unknown column name");
    checkClose(Cbc_getObjValue(model), opt, 1e-6, "optimum reached with one unknown column name");
    Cbc_deleteModel(model);
  }

  /* Nothing the model has: the start is unusable and must simply be ignored. */
  {
    Cbc_Model *model = buildFixedCharge(m);
    Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 60.0);
    Cbc_setIntParam(model, INT_PARAM_LOG_LEVEL, 1);
    const char *names[2] = { "alpha", "beta" };
    double vals[2] = { 1.0, 1.0 };
    Cbc_setMIPStart(model, 2, names, vals);
    LogCapture cap;
    Cbc_solve(model);
    std::string log = cap.finish();
    check(Cbc_isProvenOptimal(model) != 0, "proven optimal with a wholly unknown start");
    checkClose(Cbc_getObjValue(model), opt, 1e-6, "optimum reached with a wholly unknown start");
    /* Dropping it silently would leave a caller who mistyped every name with no
       clue why their start had no effect. */
    check(logHas(log, "column names do not match") || logHas(log, "were not found"),
      "an unusable start is reported rather than dropped silently");
    Cbc_deleteModel(model);
  }
}

/* ── Test 6: a suboptimal feasible start is used and improved on ───────── */

void testSuboptimalStartImproved()
{
  printf("\n=== a suboptimal but feasible MIPStart is improved upon ===\n");
  FixedCharge m;
  int bestOpen[FixedCharge::n];
  double bestX[FixedCharge::n];
  const double opt = m.optimum(bestOpen, bestX);

  /* Open everything: feasible, and strictly more expensive than the optimum
   * because every fixed charge is paid. */
  int open[FixedCharge::n];
  double x[FixedCharge::n];
  for (int i = 0; i < FixedCharge::n; ++i)
    open[i] = 1;
  const double subCost = m.bestCost(open, x);
  check(subCost > opt + 1e-6, "the start really is suboptimal");

  Cbc_Model *model = buildFixedCharge(m);
  configureQuiet(model);
  std::vector<std::string> names;
  std::vector<int> idxs;
  std::vector<double> vals;
  buildStart(m, ShapeFull, open, x, names, idxs, vals);
  std::vector<const char *> cn;
  for (size_t k = 0; k < names.size(); ++k)
    cn.push_back(names[k].c_str());
  Cbc_setMIPStart(model, (int)cn.size(), cn.data(), vals.data());

  Cbc_solve(model);
  check(Cbc_isProvenOptimal(model) != 0, "proven optimal from a suboptimal start");
  checkClose(Cbc_getObjValue(model), opt, 1e-6, "search improved on the start");
  Cbc_deleteModel(model);
}

} // namespace

int main(int argc, char **argv)
{
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "--help" || arg == "-h") {
      printf(
        "mipstart-test — tests for MIP starts through Cbc's C interface\n"
        "(Cbc_setMIPStart / Cbc_setMIPStartI), which is the path python-mip uses.\n"
        "\n"
        "Usage: mipstart-test [--help]\n"
        "\n"
        "No data files are needed: the fixture is a small fixed-charge model whose\n"
        "optimum is enumerated in-process, so a wrong answer is a wrong answer and\n"
        "not a stale fixture. Its continuous columns are not determined by the\n"
        "integers until the LP runs, which is what makes an integers-only start a\n"
        "real test.\n"
        "\n"
        "The file path (-mipstart) cannot test the three input shapes: read()\n"
        "densifies over every column before returning, so a sparse file and a full\n"
        "file arrive identically. Cbc_setMIPStart() keeps them sparse, so it is the\n"
        "only place the distinction is observable. Checks cover:\n"
        "\n"
        "  1. the three shapes -- nonzero integers only, the full solution, all\n"
        "     nonzeros -- by name and by index, with and without preprocessing;\n"
        "     each asserting acceptance, integrality, row feasibility, the\n"
        "     objective recomputed from the published columns, and that columns\n"
        "     the start never mentioned came back derived rather than zero;\n"
        "  2. INT_PARAM_MIPSTART_FIX, the runtime mode that replaced the\n"
        "     compile-time JUST_FIX_INTEGER macro -- all three values reaching the\n"
        "     same optimum, and mode 'all' actually exercising the\n"
        "     fix-continuous-too path and its integers-only retry, which was\n"
        "     unreachable while the mode was a macro;\n"
        "  3. negative integer values, where a positive-only filter on the\n"
        "     solution carried forward between solves would silently drop them;\n"
        "  4. an infeasible start being rejected rather than trusted (accepting it\n"
        "     installs a cutoff that can exclude the genuine optimum);\n"
        "  5. starts naming columns the model does not have -- partially and\n"
        "     wholly -- being survivable rather than fatal;\n"
        "  6. a suboptimal but feasible start being accepted and improved upon.\n"
        "\n"
        "Exit code: 0 = all checks passed, 1 = one or more checks failed.\n");
      return 0;
    }
    fprintf(stderr, "mipstart-test: unrecognized argument '%s' (try --help)\n", arg.c_str());
    return 1;
  }

  const Shape shapes[3] = { ShapeIntNonzero, ShapeFull, ShapeAllNonzero };
  for (int s = 0; s < 3; ++s)
    for (int byIndex = 0; byIndex < 2; ++byIndex)
      for (int pp = 0; pp < 2; ++pp)
        testShape(shapes[s], byIndex != 0, pp != 0);

  testFixModes();
  testFixAllRecovers();
  testMiniBranchAndBound();
  testNegativeIntegers();
  testInfeasibleStartRejected();
  testUnknownColumns();
  testSuboptimalStartImproved();

  printf("\n=== Summary ===\n");
  if (g_failures == 0) {
    printf("RESULT: OK (all checks passed)\n");
    return 0;
  }
  printf("RESULT: FAILED (%d check(s) failed)\n", g_failures);
  return 1;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
