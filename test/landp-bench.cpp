/**
 * Stand-alone benchmark for CglLandP.
 *
 * @file landp-bench.cpp
 * @brief replay a CglLandP call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcGomoryFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, LP solution, column types -- warm-starts the LP, and runs N
 * rounds of CglLandP. Reaching this state inside CBC costs a preprocess plus a
 * root LP, so this exists to make the iterate-measure loop on the lift-and-project
 * algorithm itself affordable: milliseconds instead of a full solve per
 * experiment.
 *
 * Figures of merit, in order of authority (BENCHMARKING-CUT-GENERATORS.md §1):
 *
 *   - **bound improvement on reoptimizing** (`objImprove`) is definitive. It is
 *     the only column that says the cuts tightened the relaxation, which is the
 *     entire point of generating them.
 *   - total *violation* (`totalViol`) is the useful proxy, available before any
 *     re-solve, and unlike a count it distinguishes fewer-but-deeper cuts.
 *   - separation *time*, which is what this exercise is trying to reduce. A gain
 *     here only counts if the bound holds.
 *   - the cut counts are reported and are the **weakest** figures: more cuts at
 *     equal bound movement is strictly worse, being more rows for the same
 *     tightening. Never rank by them and never compare them against zero.
 *
 * THE GOMORY FIXTURES ARE THE RIGHT FIXTURES, and that is a claim about LandP's
 * inputs rather than a convenience. CglLandP reads exactly four things from the
 * solver it is handed: the matrix and bounds, `isContinuous()` per column
 * (CglLandP.cpp:297), `getWarmStart()` (:236) and `getColSolution()` /
 * `getRowActivity()` (:329-330). All four are in the five-file fixture. It caches
 * nothing across calls except `numrows_`, `originalColLower_/Upper_` and
 * `extraCuts_`, all of which a fresh generator rebuilds -- so unlike Lagomory,
 * which needed a capture taken *after* cuts had been added to the LP, there is no
 * shape requirement beyond "an optimal basis of the model". Hence `gomoryFixtures/`
 * replays LandP faithfully and no separate dump was written.
 *
 * WHAT IS DIFFERENT ABOUT LandP.
 *
 * 1. THE BASIS IS THE ALGORITHM'S INPUT, NOT PROVENANCE -- more so than for
 *    Gomory. needsOptimalBasis() returns true, and LandP does not merely read a
 *    tableau row at that basis: it *pivots away from it*, up to `pivotLimit`
 *    times, choosing each pivot by a cut-improvement criterion computed from the
 *    current tableau (CglLandPSimplex.cpp:668-720). A different optimal basis of
 *    the same vertex therefore starts a different search and ends at a different
 *    cut after a different number of pivots. So `warmStartIters != 0` does not
 *    merely perturb the output, it invalidates the comparison, and the CSV
 *    carries the column so a sweep can be filtered on it.
 *
 * 2. IT CLONES THE SOLVER TWICE PER CANDIDATE ROW, and the first clone is dead.
 *    CglLandP.cpp:721 does `ncSi = t_si->clone(); landpSi.setSi(ncSi)`, and
 *    `optimize()` then opens with `delete si_; si_ = cached.solver_->clone()`
 *    (CglLandPSimplex.cpp:645-646) without having read `si_` in between --
 *    `setSi` only stores the pointer (CglLandPSimplex.hpp:78-101). So each
 *    candidate pays two full OsiClpSolverInterface deep copies where one is used.
 *    The same lines make `setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX)` and
 *    `messageHandler()->setLogLevel(0)` at :723-724 apply to an object destroyed
 *    microseconds later, so neither ever reaches the solver that actually pivots.
 *    Removing the first clone is not a one-line deletion: the clone also exists so
 *    that `optimize()`'s `delete si_` does not destroy the caller's `&si`, which
 *    the CglLandPSimplex constructor stores. `--pivot-limit=0` prices this whole
 *    path, since that branch calls `generateMig()` and clones nothing.
 *
 * 3. A RANGE ROW MAKES LandP SWITCH ITSELF OFF FOR THE REST OF THE SOLVE.
 *    On a model with any row satisfying `-1e50 < lower < upper < 1e50`, :599-668
 *    clones the solver, appends one slack column per range row, rebuilds the basis
 *    and solution and re-solves -- and then :880 executes
 *    `params_.maximumCutLength = -params_.maximumCutLength`, which the guard at
 *    :552 turns into an immediate `return` on every later call. This bench reports
 *    it exactly rather than inferring it: `rangeDisabled` is read back from
 *    `parameter().maximumCutLength < 0` after the call. Note the *other*
 *    range-row switch-off, at :548, is dead -- its condition is
 *    `if (numberRanges && false)`.
 *    Consequence for measurement: on a range-row fixture, `--rounds=2` measures
 *    one real call and one immediate return, and a fresh generator per round (which
 *    this bench does, matching CBC) hides that. Read `rangeRows` before reading a
 *    multi-round row.
 *
 * 4. CUTS CAN BE THROWN AWAY AFTER BEING GENERATED, VALIDATED AND ACCEPTED.
 *    On the range-row path a cut whose support touches an added slack column
 *    cannot be expressed in the original column space, so :810-826 skips it. That
 *    is real work producing no output, it is invisible in the cut count, and it is
 *    only reachable when `rangeRows > 0`.
 *
 * 5. THE VALIDATOR'S REJECTION COUNTERS ARE A DIAGNOSTIC, NOT A CENSUS.
 *    `validator().numRejected(code)` is reported per code because it says *why*
 *    the pivoting produced nothing, which no other column does. But
 *    CglLandPValidator.cpp:56-57 returns `SmallViolation` WITHOUT incrementing the
 *    counter (the increment at :167 is a second, later path to the same code), so
 *    `rejSmallViol` undercounts. Do not compute an acceptance rate from these.
 *
 * 6. `info.pass` IS NOT A BUDGET SWITCH HERE. Unlike Gomory, LandP reads only
 *    `info.inTree` (:562-566, lowering pivotLimit to pivotLimitInTree) and
 *    `info.pass` purely for a log message (:558). So `--pass=N` is cosmetic and
 *    `--in-tree` is the real mode switch. The default 0/false is CBC's root call.
 *
 * A THING THAT LOOKS LIKE A DEFECT AND IS NOT, recorded so it is not chased
 * twice. `params_.timeLimit += CoinCpuTime()` at :705 and `-= CoinCpuTime()` at
 * :867 look like a botched conversion to an absolute deadline, because the local
 * `params` copy taken at :556 keeps the relative value. They are not: `optimize()`
 * does its own `timeLimit = min(params.timeLimit, params.singleCutTimeLimit) +
 * CoinCpuTime()` (CglLandPSimplex.cpp:627-628), so the relative value is what it
 * wants. The pair on the member is a separate mechanism -- it debits the member by
 * the CPU the round consumed, so `params.timeLimit < 0` at :565 switches pivoting
 * off once the cross-call budget is spent. Inert at the COIN_DBL_MAX default.
 *
 * THE CONTROL FLAGS. §7 of the process doc: price a stage before optimizing it.
 * LandP admits a three-point decomposition, which no other generator here does:
 *
 *   --away=0.5           the floor. `away` gates candidate selection at :918-931
 *                        (`INT_INFEAS(x) <= away` is skipped), so 0.5 admits only
 *                        exactly-half-integral columns. What remains timed is the
 *                        irreducible per-call cost: `cached_.getData()` including
 *                        its one solver clone (:355), the CglLandPSimplex
 *                        constructor, and the candidate sort. time(0.5) is the
 *                        floor no change to the candidate loop can beat.
 *                        CHECK rowCuts_n BEFORE BELIEVING IT -- the test is
 *                        strict, so a basic integer at exactly 0.5 still passes at
 *                        away=0.5, and half-integral LPs (graph colouring) admit
 *                        nearly everything.
 *   --pivot-limit=0      the no-pivoting path: `generateMig()` per candidate, no
 *                        clones, no cut-improving simplex. So
 *                        time(default) - time(pl=0) prices the entire L&P search
 *                        including both per-candidate clones, and
 *                        time(pl=0) - time(away=0.5) prices the per-candidate MIG
 *                        construction alone. THE CUTS DIFFER at pivotLimit=0 --
 *                        they are plain mixed-integer Gomory cuts -- so this is a
 *                        timing control only, never an exactness baseline.
 *
 * Usage:
 *   landp-bench <stem> [options]           (stem = dir/name.tag)
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.bas, <stem>.sol,
 * <stem>.ctype, <stem>.meta. A full path to any one of those also works; the
 * suffix is stripped.
 */

#include "CglLandP.hpp"
#include "CglLandPValidator.hpp"
#include "CglTreeInfo.hpp"
#include "ClpSimplex.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>

/// Wall-clock, not CPU: separation cost as a user would feel it.
static double wallClock()
{
  return CoinGetTimeOfDay();
}

/// CglLandP.cpp:14. Duplicated rather than included because the macro is defined
/// in the .cpp, not the header; benchCandidates below must use the same rule.
#define LANDP_INT_INFEAS(value) fabs((value) - floor((value) + 0.5))

/*
 * Per-stage attribution, only when Cgl was built -DCGL_LANDP_PROFILE.
 *
 * Declared here rather than in CglLandP.hpp on purpose: this is a local
 * diagnostic build, not API, and putting it in the shipped header would mean a
 * header/library mismatch produces a link error for anyone who builds Cgl
 * normally. The cost of that choice is that these two prototypes must match the
 * definitions in CglLandP.cpp by hand.
 *
 * The accumulators are file-static, so a profile run must be serial -- which is
 * already the rule for any run whose times get quoted.
 */
#ifdef CGL_LANDP_PROFILE
void cglLandPProfileReset();
void cglLandPProfilePrint(const char *tag);
#define LANDP_PROF_RESET() cglLandPProfileReset()
#define LANDP_PROF_PRINT(tag) cglLandPProfilePrint(tag)
#else
#define LANDP_PROF_RESET()
#define LANDP_PROF_PRINT(tag)
#endif

static bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

/**
 * Reduce any of the fixture's file names to the shared stem, so a caller can pass
 * whichever one tab-completion produced.
 */
static std::string fixtureStem(const char *arg)
{
  std::string s(arg);
  static const char *suffixes[]
    = { ".mps.gz", ".mps", ".bas", ".sol", ".ctype", ".meta", ".bas.status" };
  for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
    const std::string suf(suffixes[i]);
    if (s.size() > suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0)
      return s.substr(0, s.size() - suf.size());
  }
  return s;
}

/// Last path component, for the CSV name column.
static std::string baseName(const std::string &path)
{
  const size_t slash = path.rfind('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

/**
 * Read one numeric key out of the fixture's `.meta`. Returns `dflt` when the file
 * or the key is absent, so a fixture written before `.meta` carried that key
 * still loads.
 */
static double metaNum(const std::string &path, const char *key, double dflt)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return dflt;

  char line[512];
  double value = dflt;
  while (fgets(line, sizeof(line), fp)) {
    char k[256];
    double v = 0.0;
    if (sscanf(line, "%255s %lf", k, &v) == 2 && strcmp(k, key) == 0) {
      value = v;
      break;
    }
  }
  fclose(fp);
  return value;
}

static long metaInt(const std::string &path, const char *key, long dflt)
{
  return (long)metaNum(path, key, (double)dflt);
}

/**
 * Remove the pad row, if this fixture has one.
 *
 * A capture whose matrix held empty columns carries one extra final row, without
 * which writeMps would drop those columns and shift every index after them -- see
 * cbcGomoryFixtureWriteMps. For LandP the row is not merely cosmetic even though
 * it is redundant: it is an extra basis position, so leaving it in changes both
 * the tableau LandP pivots on and the `nBasics_` cap on the cut count, and it
 * would also make the `.bas` one artificial short of the matrix and so unusable.
 *
 * Both conditions are required before deleting anything: `.meta` must say padding
 * happened, *and* the loaded row count must be exactly one more than the captured
 * one. Either alone could delete a real row from a fixture whose meta is stale.
 */
static bool dropPadRow(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string meta = stem + ".meta";
  if (metaInt(meta, "paddedColumns", 0) <= 0)
    return false;

  const long capturedRows = metaInt(meta, "rows", -1);
  if (capturedRows < 0 || si.getNumRows() != (int)capturedRows + 1) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: meta says padded but rows=%d against captured %ld; "
                      "leaving the matrix alone\n",
        baseName(stem).c_str(), si.getNumRows(), capturedRows);
    return false;
  }

  const int last = si.getNumRows() - 1;
  si.deleteRows(1, &last);
  return true;
}

/**
 * Restore integrality from the `.ctype` sidecar.
 *
 * MPS conveys integrality only through the bound type, and a column with
 * lb == ub takes writeMps's " FX " branch, which has no integer form -- so every
 * integer column CBC had fixed by bound tightening reads back **continuous**.
 *
 * For LandP a lost marker changes the experiment in three places. `integers_[]`
 * is built from `isContinuous()` at CglLandP.cpp:297 and then gates candidate
 * selection at :922, so a lost marker removes the column from the search
 * entirely. It also propagates: :308-326 clears `integerSlacks[]` for every row
 * a *continuous* column touches, so one lost marker can disqualify whole rows'
 * slacks from contributing integral coefficients. And CglLandPSimplex reads the
 * same array when it builds the cut coefficient for each nonbasic
 * (`integers_[]`, the mixed-integer versus continuous formula), so the cuts that
 * *are* generated move too. Skipping this step does not merely measure a smaller
 * problem, it measures different cuts.
 *
 * Returns the number of columns re-marked, or -1 when no usable sidecar was
 * found. A sidecar whose column count disagrees with the model is refused rather
 * than partly applied: marking arbitrary columns integer is worse than the loss
 * it repairs.
 */
static int restoreColTypes(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string path = stem + ".ctype";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: no .ctype sidecar; integer columns that were "
                      "fixed at capture will read back continuous\n",
        baseName(stem).c_str());
    return -1;
  }

  int sidecarCols = -1;
  if (fscanf(fp, "cols %d\n", &sidecarCols) != 1 || sidecarCols != si.getNumCols()) {
    fprintf(stderr, "ERROR: %s: .ctype is for %d columns, model has %d; ignoring it\n",
      baseName(stem).c_str(), sidecarCols, si.getNumCols());
    fclose(fp);
    return -1;
  }

  int idx = 0, type = 0, restored = 0;
  while (fscanf(fp, "%d %d\n", &idx, &type) == 2) {
    if (idx < 0 || idx >= si.getNumCols()) {
      fprintf(stderr, "ERROR: %s: .ctype names column %d, out of range\n",
        baseName(stem).c_str(), idx);
      fclose(fp);
      return -1;
    }
    if (si.isContinuous(idx)) {
      si.setInteger(idx);
      ++restored;
    }
  }
  fclose(fp);

  // getColType() caches. CglLandP reads integrality through isContinuous(), which
  // OsiClp answers from the same cache, so it must be recomputed or the generator
  // sees the pre-restore view.
  si.getColType(true);
  return restored;
}

/*
 * Emit every cut, coefficient by coefficient, in a form where textual equality
 * is bit equality.
 *
 * The CSV fields are a fingerprint, not a proof: totalViol and avgCutLen are sums
 * over all cuts, so two different cut sets can agree on every one of them. When
 * the claim being checked is "this optimization changes nothing", a fingerprint
 * match is the weaker statement. %a prints the exact IEEE double, so a one-ulp
 * difference in a single coefficient of a single cut shows up as a diff rather
 * than rounding away at the 15th digit.
 *
 * Order is as generated, deliberately unsorted: the order cuts land in OsiCuts is
 * itself part of what must not change, since applyCuts consumes them in that
 * order and later rounds depend on it. For LandP there is a second reason -- the
 * candidate order is `getSortedFractionalIndices`' sort, and a change that
 * reorders it produces the same cuts in a different order, which is a behaviour
 * change worth seeing rather than normalizing away.
 */
static void dumpCuts(const OsiCuts &cs, const std::string &name, int round)
{
  for (int c = 0; c < cs.sizeRowCuts(); ++c) {
    const OsiRowCut &rc = cs.rowCut(c);
    const CoinPackedVector &row = rc.row();
    const int len = row.getNumElements();
    const int *idx = row.getIndices();
    const double *el = row.getElements();
    printf("[cut] %s r%d row %d n %d lb %a ub %a ge %a", name.c_str(), round, c,
      len, rc.lb(), rc.ub(), rc.globallyValid() ? 1.0 : 0.0);
    for (int k = 0; k < len; ++k)
      printf(" %d:%a", idx[k], el[k]);
    printf("\n");
  }
  // LandP emits no column cuts -- every insert in generateCuts is
  // insertIfNotDuplicate of an OsiRowCut -- so this loop should never print.
  // Kept so that a change which starts producing them is visible rather than
  // silently dropped from the exactness gate.
  for (int c = 0; c < cs.sizeColCuts(); ++c) {
    const OsiColCut &cc = cs.colCut(c);
    const CoinPackedVector &lbs = cc.lbs();
    const CoinPackedVector &ubs = cc.ubs();
    printf("[cut] %s r%d col %d nlb %d nub %d", name.c_str(), round, c,
      lbs.getNumElements(), ubs.getNumElements());
    for (int k = 0; k < lbs.getNumElements(); ++k)
      printf(" L%d:%a", lbs.getIndices()[k], lbs.getElements()[k]);
    for (int k = 0; k < ubs.getNumElements(); ++k)
      printf(" U%d:%a", ubs.getIndices()[k], ubs.getElements()[k]);
    printf("\n");
  }
}

/// Everything the fixture determines, plus how long loading it took.
struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  /// Pivots the warm start needed. Zero is the expected value, and for LandP it
  /// is not merely a hygiene check -- see loadFixture.
  int warmStartIters = 0;
  bool paddedRowDropped = false;
  int restoredColTypes = -1;
  bool haveBasis = false;
  bool ok = false;
};

/**
 * Load the fixture and warm-start to the captured optimum.
 *
 * Two things are needed to actually *land* on the captured vertex, and getting
 * either wrong looks like success while silently changing the experiment. For
 * LandP the stakes are the highest of any generator benched here: the basis is
 * not just the point the cut is read off, it is the starting point of a search
 * that pivots away from it, so landing on a different optimal basis of the same
 * vertex sends the search down a different path and produces a different,
 * equally valid cut after a different amount of work. A before/after comparison
 * across such a load is meaningless, not merely noisy.
 *
 * First, `readBasis` writes into ClpSimplex's `status_`, but OsiClp caches a
 * separate `CoinWarmStartBasis basis_` and `resolve()` overwrites the model from
 * it. `setWarmStart(NULL)` refreshes that cache from the model, which is what
 * makes the file's basis survive into the solve. It is also what makes
 * `si.getWarmStart()` -- which is where CglLandP reads the basis, at
 * CglLandP.cpp:236 -- return the captured one.
 *
 * Second, `resolve()` rather than `initialSolve()`: presolve discards the basis.
 * With both in place an already-optimal fixture costs 0 iterations, which is the
 * cheap self-check that the warm start worked. `setPerturbation(50)` matches the
 * generator; perturbation left on moves the vertex even from a correct basis.
 *
 * A fixture with no usable basis is refused outright rather than solved cold.
 */
static bool loadFixture(Fixture &f, const std::string &stem, bool quiet)
{
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";

  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no problem file for stem %s\n", stem.c_str());
    return false;
  }

  ClpSimplex *lp = f.si.getModelPtr();
  lp->setLogLevel(0);
  if (f.si.readMps(mps.c_str())) {
    fprintf(stderr, "ERROR: failed to read %s\n", mps.c_str());
    return false;
  }

  // Before the basis: the basis was written from the captured model, so it has
  // one artificial fewer than the padded matrix and would not line up.
  f.paddedRowDropped = dropPadRow(f.si, stem, quiet);

  // Before the generator: integrality gates candidate selection, the integer-slack
  // classification and the cut coefficient formula.
  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  if (!fileExists(bas)) {
    fprintf(stderr, "ERROR: %s: no .bas. LandP pivots away from a specific optimal "
                    "basis, so there is nothing to replay without it\n",
      baseName(stem).c_str());
    return false;
  }
  if (lp->readBasis(bas.c_str()) < 0) {
    fprintf(stderr, "ERROR: %s: failed to read basis %s\n",
      baseName(stem).c_str(), bas.c_str());
    return false;
  }
  f.haveBasis = true;
  // Push the freshly-read status into OsiClp's cached basis_, or the solve below
  // installs the stale cache over it -- and getWarmStart() would hand CglLandP
  // that stale cache too.
  f.si.setWarmStart(NULL);

  lp->setPerturbation(50);

  const double t0 = wallClock();
  f.si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
  f.si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
  f.si.resolve();
  f.warmStartTime = wallClock() - t0;
  f.warmStartIters = f.si.getIterationCount();

  if (!f.si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: LP not optimal after warm start (%s)\n", stem.c_str());
    return false;
  }
  // A correct warm start from an optimal basis costs no pivots. Anything else
  // means the fixture landed on a different vertex or a different basis of the
  // same one, which for LandP is a different search -- so the CSV carries
  // warmStartIters and a sweep should filter on it.
  if (f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured basis "
                    "did not survive, so these are cuts from a DIFFERENT basis and "
                    "not comparable across runs\n",
      baseName(stem).c_str(), f.warmStartIters);
  }

  f.ok = true;
  return true;
}

/**
 * Count the rows LandP would treat as ranges, by CglLandP.cpp:541-546's rule.
 *
 * Reported because a nonzero value changes the call in four ways at once: the
 * solver is cloned and augmented with one slack column per range row, the basis
 * and solution are rebuilt, accepted cuts touching an added slack are silently
 * discarded, and the generator switches itself off for every subsequent call.
 * Any row of the CSV with rangeRows > 0 is measuring a different code path.
 */
static int countRangeRows(const OsiSolverInterface &si)
{
  const double *lo = si.getRowLower();
  const double *up = si.getRowUpper();
  int n = 0;
  for (int i = 0; i < si.getNumRows(); ++i)
    if (lo[i] < up[i] && lo[i] > -1.0e50 && up[i] < 1.0e50)
      ++n;
  return n;
}

/**
 * How many candidate rows LandP will iterate over, by the rule at
 * CglLandP.cpp:918-931 (basic, structural, integer, INT_INFEAS > away), capped by
 * `maximumCandidates`.
 *
 * Replicated here rather than read from the generator because the count is not
 * exposed. It is the denominator for every per-candidate cost, and it is what
 * makes `--away=0.5` legible: the flag is only a control if this number actually
 * collapses, which on half-integral LPs it does not.
 *
 * The rule is duplicated, so it can drift from the generator's. It is
 * informational for that reason -- never gate a conclusion on it without checking
 * rowCuts_n moved the same way.
 */
static int benchCandidates(const OsiSolverInterface &si, double away, int cap)
{
  const CoinWarmStartBasis *ws
    = dynamic_cast< const CoinWarmStartBasis * >(si.getWarmStart());
  if (!ws)
    return -1;
  const double *x = si.getColSolution();
  const int n = si.getNumCols();
  int count = 0;
  for (int j = 0; j < n && j < ws->getNumStructural(); ++j) {
    if (ws->getStructStatus(j) != CoinWarmStartBasis::basic)
      continue;
    if (si.isContinuous(j))
      continue;
    if (LANDP_INT_INFEAS(x[j]) <= away)
      continue;
    ++count;
  }
  delete ws;
  return count < cap ? count : cap;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.bas/.sol/.ctype/.meta are appended.\n"
    "Passing any one of those files also works. The .bas is REQUIRED: LandP pivots\n"
    "away from a specific optimal basis. gomoryFixtures/ replays LandP faithfully --\n"
    "see the file comment for why no separate dump was needed.\n"
    "\n"
    "Call shape (defaults reproduce CBC's root call -- see the file comment):\n"
    "  --rounds=N          LandP rounds, LP re-solved between them (default 1;\n"
    "                      CBC calls it once per pass with a fresh generator)\n"
    "  --pass=N            info.pass (default 0). COSMETIC for LandP: the only\n"
    "                      reader is a log message at CglLandP.cpp:558. Unlike\n"
    "                      Gomory, pass does not change any budget or tolerance.\n"
    "  --in-tree           info.inTree = true, the real mode switch: :562-566\n"
    "                      lowers pivotLimit to pivotLimitInTree and turns on\n"
    "                      countMistakenRc.\n"
    "  --options=N         info.options bitmask (default 0). LandP reads it only\n"
    "                      through CglCutGenerator, so it is here for symmetry.\n"
    "\n"
    "CBC's knobs (defaults are CBC's, from CbcSolverCutSetup.cpp:431-450, NOT\n"
    "CglLandP's constructor defaults, which differ on two of them):\n"
    "  --maximum-cut-length=N  cut length cap (default 2000 = CBC's; CglLandP's own\n"
    "                      default is 10000, and -lift onglobal uses 2000000).\n"
    "                      NOTE a NEGATIVE value here is not a smaller cap -- :552\n"
    "                      reads negative as \"switched off\" and returns at once.\n"
    "  --min-violation=F   validator's minimum violation (default 1e-4 = CBC's;\n"
    "                      CglLandP's own default is 0, i.e. accept any violation)\n"
    "\n"
    "LandP's own knobs (defaults from CglLandP::Parameters::Parameters()):\n"
    "  --pivot-limit=N     cut-improving pivots per candidate (default 20).\n"
    "                      0 IS A CONTROL FLAG, not a tighter budget: it takes the\n"
    "                      generateMig() branch, which clones nothing and does no\n"
    "                      pivoting, so time(default)-time(0) prices the whole L&P\n"
    "                      search including both per-candidate solver clones. The\n"
    "                      CUTS DIFFER at 0 (plain MIG), so never use it as an\n"
    "                      exactness baseline.\n"
    "  --pivot-limit-in-tree=N  (default 10) only bites with --in-tree\n"
    "  --away=F            fractionality gate (default 5e-4). Usable as a control\n"
    "                      flag at 0.5, where almost no candidate is admitted, so\n"
    "                      time(0.5) approximates the irreducible per-call floor\n"
    "                      (getData + its clone, the simplex-wrapper constructor,\n"
    "                      the candidate sort). CHECK rowCuts_n AND benchCandidates\n"
    "                      BEFORE BELIEVING IT: the test is strict, so a basic\n"
    "                      integer sitting at exactly 0.5 passes at any away<=0.5.\n"
    "  --max-cut-per-round=N   (default 5000). The loop also stops at nBasics_\n"
    "                      (= number of rows), whichever binds first.\n"
    "  --maximum-candidates=N  (default 1000000) truncates the sorted candidate\n"
    "                      list, so it caps work directly.\n"
    "  --pivot-tol=F       (default 1e-4)\n"
    "  --failed-pivot-limit=N      (default 1)\n"
    "  --degenerate-pivot-limit=N  (default 0)\n"
    "  --extra-cuts-limit=N        (default 5)\n"
    "  --extra-cuts=N      0 none (default), 1 AtOptimalBasis, 2 WhenEnteringBasis,\n"
    "                      3 AllViolatedMigs. 3 also runs genThisBasisMigs per\n"
    "                      candidate, so it is a large behaviour change.\n"
    "  --pivot-selection=N 0 mostNegativeRc (default), 1 bestPivot,\n"
    "                      2 initialReducedCosts\n"
    "  --sep-space=N       0 Fractional (default), 1 Fractional_rc, 2 Full\n"
    "  --normalization=N   0 Unweighted (default), 1 WeightRHS, 2 WeightLHS,\n"
    "                      3 WeightBoth\n"
    "  --lhs-norm=N        0 L1 (default), 1 L2, 2 SupportSize, 3 Infinity,\n"
    "                      4 Average, 5 Uniform\n"
    "  --rhs-weight-type=N 0 Fixed (default), 1 Dynamic. NOTE generateCuts\n"
    "                      overwrites rhsWeight with numrows+2 at :557 regardless.\n"
    "  --modularize        Balas's modularization (default off)\n"
    "  --no-strengthen     turn off final-cut strengthening (default on)\n"
    "  --no-perturb        turn off the perturbation procedure (default on)\n"
    "  --count-mistaken-rc charge a mistaken reduced cost against failedPivotLimit\n"
    "                      (default off at the root, and --in-tree forces it ON\n"
    "                      together with pivotLimitInTree). This flag exists to\n"
    "                      separate those two: with perturbation on, ~96%% of row\n"
    "                      candidates are mistaken, and at the root nothing bounds\n"
    "                      the rescan loop that chases them. NOT output-neutral --\n"
    "                      it can cost pivots, so read objImprove, not just time.\n"
    "  --no-exact-retry    pick a row's pivot candidate from the tabulated reduced\n"
    "                      costs, as upstream does, instead of from the exact cost\n"
    "                      of all four of its (direction, gammaSign) candidates\n"
    "                      computed from the tableau row already in hand (default:\n"
    "                      exact). This is the control arm: with it, counters must\n"
    "                      match a build from before the exact route existed,\n"
    "                      exactly. Without --exact-best the exact cost is only a\n"
    "                      screen on the candidate the tables named, so the pivots\n"
    "                      and the cuts are unchanged and only the futile column\n"
    "                      searches go away.\n"
    "  --exact-best        let the exact cost pick which of a row's four candidates\n"
    "                      to pivot on rather than only confirm the tables' pick\n"
    "                      (default off). NOT output-neutral: it changes the pivot\n"
    "                      sequence, so read objImprove and totalViol, not time.\n"
    "  --no-pre-length-gate  clone the LP before testing the maximumCutLength\n"
    "                      gate, as upstream does, instead of reading the source\n"
    "                      row from the cached optimal-basis solver first\n"
    "                      (default: gate first). 71.7%% of optimize() calls stop\n"
    "                      at that gate having pivoted zero times, so most of\n"
    "                      the clones are built and destroyed for nothing. The\n"
    "                      cut a passing row yields is unchanged -- this is the\n"
    "                      control arm for proving that.\n"
    "  --no-tableau-row    use the disjunction rather than the tableau row\n"
    "  --time-limit=F      seconds per round (default COIN_DBL_MAX). NOT a wall to\n"
    "                      measure against: exceeding it makes later candidates take\n"
    "                      the no-pivot path, so it changes the cuts.\n"
    "  --single-cut-time-limit=F   seconds per candidate (default COIN_DBL_MAX)\n"
    "  --log-level=N       generator log level: 1 begin/end, 2 per cut. Useful for\n"
    "                      diagnosis; it writes to stdout, so not with --dump-cuts.\n"
    "\n"
    "Output:\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --self-test         run internal consistency checks and exit\n"
    "  --dump-cuts         print every cut coefficient as an exact IEEE hex\n"
    "                      double ([cut] lines on stdout, before the CSV row).\n"
    "                      Use this, not the CSV fields, to check that a change\n"
    "                      really is output-neutral: the CSV aggregates and two\n"
    "                      different cut sets can share every aggregate.\n"
    "  --no-bound          skip the post-separation LP re-solve. TIMING AID ONLY:\n"
    "                      it makes objEnd/objImprove/boundMoved meaningless (they\n"
    "                      report no movement), so never judge a change with it --\n"
    "                      the bound is the metric that decides. It exists because\n"
    "                      the re-solve can dwarf what is being measured.\n"
    "                      Refused with --rounds>1, where round 2 would otherwise\n"
    "                      separate against an LP the round-1 cuts never entered.\n"
    "  --quiet             suppress warnings\n",
    prog);
}

/**
 * Did the bound actually move, or is this floating-point noise?
 *
 * objImprove is a difference of two LP objectives, so it inherits their absolute
 * error, which scales with their magnitude -- an LP around 1e6 carries roughly
 * 1e-10 of slop. Reading such a difference against zero therefore reports an
 * improvement on nearly every fixture. A relative test separates real movement
 * from noise, with an absolute floor so an objective near zero does not make the
 * relative test hypersensitive.
 */
static bool boundMoved(double objStart, double objImprove)
{
  const double scale = fabs(objStart) > 1.0 ? fabs(objStart) : 1.0;
  return objImprove > 1.0e-9 * scale && objImprove > 1.0e-9;
}

/**
 * Internal checks that do not need a fixture, so a build can be validated before
 * any measurement. Each one guards a mistake that would otherwise surface as a
 * plausible-looking number.
 */
static int selfTest()
{
  int failures = 0;

  // fixtureStem must strip every suffix the dump writes, and must not eat a
  // directory name that happens to end in one of them. The tag is .gomory here
  // deliberately: these ARE the Gomory fixtures.
  struct {
    const char *in;
    const char *want;
  } stems[] = {
    { "d/x.gomory.mps.gz", "d/x.gomory" },
    { "d/x.gomory.bas", "d/x.gomory" },
    { "d/x.gomory.bas.status", "d/x.gomory" },
    { "d/x.gomory.meta", "d/x.gomory" },
    { "d/x.gomory", "d/x.gomory" },
  };
  for (size_t i = 0; i < sizeof(stems) / sizeof(stems[0]); ++i) {
    const std::string got = fixtureStem(stems[i].in);
    if (got != stems[i].want) {
      fprintf(stderr, "FAIL fixtureStem(%s) = %s, want %s\n",
        stems[i].in, got.c_str(), stems[i].want);
      ++failures;
    }
  }

  // boundMoved must be scale-relative: 1e-10 on an objective of 1e6 is noise,
  // the same 1e-10 on an objective of 1 is not necessarily, and a real move on a
  // large objective must still register.
  if (boundMoved(1.0e6, 1.0e-10)) {
    fprintf(stderr, "FAIL boundMoved: 1e-10 on 1e6 should be noise\n");
    ++failures;
  }
  if (!boundMoved(1.0e6, 1.0e-2)) {
    fprintf(stderr, "FAIL boundMoved: 1e-2 on 1e6 should count\n");
    ++failures;
  }
  if (boundMoved(0.0, 1.0e-12)) {
    fprintf(stderr, "FAIL boundMoved: 1e-12 on 0 should be noise\n");
    ++failures;
  }

  // LANDP_INT_INFEAS is a copy of the generator's macro, and benchCandidates
  // depends on it agreeing. Pin the two properties that matter: it is distance to
  // the NEAREST integer (not the floor), and it is symmetric.
  if (fabs(LANDP_INT_INFEAS(2.9) - 0.1) > 1.0e-15
    || fabs(LANDP_INT_INFEAS(3.1) - 0.1) > 1.0e-15
    || fabs(LANDP_INT_INFEAS(3.0)) > 1.0e-15) {
    fprintf(stderr, "FAIL LANDP_INT_INFEAS is not distance-to-nearest-integer: "
                    "%g %g %g\n",
      LANDP_INT_INFEAS(2.9), LANDP_INT_INFEAS(3.1), LANDP_INT_INFEAS(3.0));
    ++failures;
  }

  // The knobs this bench sets are PLAIN PUBLIC MEMBERS of CglLandP::Parameters,
  // not setters with validation -- so unlike CglGomory there is no silent-clamp
  // trap. Pin that: if any of them ever grows a setter that rejects a value, the
  // CSV would report what was asked for while the run used something else.
  {
    CglLandP g;
    g.parameter().pivotLimit = 0;
    g.parameter().away = 0.5;
    g.parameter().maximumCutLength = 2000;
    g.parameter().maxCutPerRound = 5000;
    g.validator().setMinViolation(1.0e-4);
    if (g.parameter().pivotLimit != 0 || fabs(g.parameter().away - 0.5) > 1.0e-15
      || g.parameter().maximumCutLength != 2000
      || g.parameter().maxCutPerRound != 5000
      || fabs(g.validator().getMinViolation() - 1.0e-4) > 1.0e-15) {
      fprintf(stderr, "FAIL the knobs do not round-trip: pivotLimit=%d away=%g "
                      "maximumCutLength=%d maxCutPerRound=%d minViolation=%g\n",
        g.parameter().pivotLimit, g.parameter().away,
        g.parameter().maximumCutLength, g.parameter().maxCutPerRound,
        g.validator().getMinViolation());
      ++failures;
    }
  }

  // CBC's two overrides must differ from CglLandP's own defaults, or the comment
  // claiming this bench reproduces CBC's root call rather than the library default
  // is vacuous.
  {
    CglLandP g;
    if (g.parameter().maximumCutLength == 2000) {
      fprintf(stderr, "FAIL CglLandP's default maximumCutLength is now 2000, the same "
                      "as CBC's override; the --maximum-cut-length docs are stale\n");
      ++failures;
    }
    if (fabs(g.validator().getMinViolation() - 1.0e-4) < 1.0e-15) {
      fprintf(stderr, "FAIL CglLandP's default minViolation is now 1e-4, the same as "
                      "CBC's override; the --min-violation docs are stale\n");
      ++failures;
    }
  }

  // needsOptimalBasis() is the premise of this whole bench -- it is why the .bas
  // is required rather than optional, and why a nonzero warmStartIters
  // invalidates a comparison. Pin it.
  {
    CglLandP g;
    if (!g.needsOptimalBasis()) {
      fprintf(stderr, "FAIL needsOptimalBasis() is now false; this bench's insistence "
                      "on the captured basis needs rethinking\n");
      ++failures;
    }
  }

  // The range-row self-disable is how rangeDisabled is read back, and the file
  // comment tells the reader to trust that column. Pin the guard's shape: a
  // negative maximumCutLength must make generateCuts return early. Checked in the
  // source rather than by running, because triggering it needs a range-row model
  // and two calls.
  {
    const char *src = "../../Cgl/src/CglLandP/CglLandP.cpp";
    FILE *fp = fopen(src, "r");
    if (!fp) {
      printf("  (skipped range-disable check: %s not readable from cwd)\n", src);
    } else {
      char line[512];
      bool guard = false, flip = false;
      // The guard spans two lines -- `} else if (params_.maximumCutLength < 0) {`
      // then `return;` -- so match the test and then look ahead a couple of lines
      // for the return rather than demanding both on one line.
      int lookahead = 0;
      while (fgets(line, sizeof(line), fp)) {
        if (lookahead > 0) {
          if (strstr(line, "return"))
            guard = true;
          --lookahead;
        }
        if (strstr(line, "params_.maximumCutLength < 0"))
          lookahead = 2;
        if (strstr(line, "params_.maximumCutLength = -params_.maximumCutLength"))
          flip = true;
      }
      fclose(fp);
      if (!guard) {
        fprintf(stderr, "FAIL no 'params_.maximumCutLength < 0 ... return' guard in "
                        "%s; rangeDisabled no longer means what the docs say\n",
          src);
        ++failures;
      }
      if (!flip) {
        fprintf(stderr, "FAIL no 'params_.maximumCutLength = -params_.maximumCutLength' "
                        "in %s; the range-row self-disable is gone and rangeDisabled "
                        "will always read 0\n",
          src);
        ++failures;
      }
    }
  }

  printf("self-test: %d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}

static const char *CSV_HEADER
  = "name,pivotLimit,pivotLimitInTree,away,maxCutPerRound,maximumCandidates,"
    "maximumCutLength,minViolation,pivotTol,failedPivotLimit,degeneratePivotLimit,"
    "extraCutsLimit,extraCuts,pivotSelection,sepSpace,normalization,lhsNorm,"
    "rhsWeightType,modularize,strengthen,perturb,useTableauRow,timeLimit,"
    "singleCutTimeLimit,pass,inTree,options,rounds,"
    "rowCuts_n,colCuts_n,rowsAdded,totalViol,maxViol,avgCutLen,maxCutLen,"
    "sepTime,warmStartTime,warmStartIters,resolveTime,resolveIters,"
    "objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "rows,cols,intCols,elements,rangeRows,benchCandidates,"
    "rejSmallViol,rejSmallCoef,rejBigDynamic,rejDense,rejEmpty,rangeDisabled,"
    "restoredInt,padDropped,rowCutsPerRound,objImprovePerRound";

int main(int argc, char *argv[])
{
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }
  if (strcmp(argv[1], "--header") == 0) {
    printf("%s\n", CSV_HEADER);
    return 0;
  }
  if (strcmp(argv[1], "--self-test") == 0)
    return selfTest();
  if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
    usage(argv[0]);
    return 0;
  }

  bool csvHeader = false;
  bool quiet = false;
  bool dumpCutsFlag = false;
  bool noBound = false;
  int logLevel = 0;
  // One round by default: CBC constructs a fresh CglLandP per cut pass, so a
  // single round is the call being optimized.
  int maxRounds = 1;

  // Start from the library defaults and overwrite only what CBC overwrites, so
  // this stays correct if CglLandP::Parameters gains a member.
  CglLandP::Parameters defaults;
  int pivotLimit = defaults.pivotLimit;
  int pivotLimitInTree = defaults.pivotLimitInTree;
  double away = defaults.away;
  int maxCutPerRound = defaults.maxCutPerRound;
  int maximumCandidates = defaults.maximumCandidates;
  double pivotTol = defaults.pivotTol;
  int failedPivotLimit = defaults.failedPivotLimit;
  int degeneratePivotLimit = defaults.degeneratePivotLimit;
  int extraCutsLimit = defaults.extraCutsLimit;
  int extraCuts = (int)defaults.generateExtraCuts;
  int pivotSelection = (int)defaults.pivotSelection;
  int sepSpace = (int)defaults.sepSpace;
  int normalization = (int)defaults.normalization;
  int lhsNorm = (int)defaults.lhs_norm;
  int rhsWeightType = (int)defaults.rhsWeightType;
  bool modularize = defaults.modularize;
  bool strengthen = defaults.strengthen;
  bool perturb = defaults.perturb;
  /* Not folded into --in-tree: that switch moves pivotLimit as well, so a run with
     it cannot say which of the two changes mattered. */
  bool countMistakenRc = defaults.countMistakenRc;
  bool exactRetry = defaults.exactRetry;
  bool exactBest = defaults.exactBest;
  bool preLengthGate = defaults.preLengthGate;
  bool useTableauRow = defaults.useTableauRow;
  double timeLimit = defaults.timeLimit;
  double singleCutTimeLimit = defaults.singleCutTimeLimit;
  // CBC's two overrides (CbcSolverCutSetup.cpp:436-437), which differ from the
  // library defaults -- selfTest pins that they still differ.
  int maximumCutLength = 2000;
  double minViolation = 1.0e-4;

  int pass = 0;
  bool inTree = false;
  int options = 0;

  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (strcmp(a, "--csv-header") == 0) {
      csvHeader = true;
    } else if (strcmp(a, "--dump-cuts") == 0) {
      dumpCutsFlag = true;
    } else if (strcmp(a, "--quiet") == 0) {
      quiet = true;
    } else if (strcmp(a, "--no-bound") == 0) {
      noBound = true;
    } else if (strcmp(a, "--in-tree") == 0) {
      inTree = true;
    } else if (strcmp(a, "--modularize") == 0) {
      modularize = true;
    } else if (strcmp(a, "--no-strengthen") == 0) {
      strengthen = false;
    } else if (strcmp(a, "--no-perturb") == 0) {
      perturb = false;
    } else if (strcmp(a, "--count-mistaken-rc") == 0) {
      countMistakenRc = true;
    } else if (strcmp(a, "--no-exact-retry") == 0) {
      exactRetry = false;
    } else if (strcmp(a, "--exact-best") == 0) {
      exactBest = true;
    } else if (strcmp(a, "--no-pre-length-gate") == 0) {
      preLengthGate = false;
    } else if (strcmp(a, "--no-tableau-row") == 0) {
      useTableauRow = false;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--pass=", 7) == 0) {
      pass = atoi(a + 7);
    } else if (strncmp(a, "--options=", 10) == 0) {
      options = atoi(a + 10);
    } else if (strncmp(a, "--log-level=", 12) == 0) {
      logLevel = atoi(a + 12);
    } else if (strncmp(a, "--pivot-limit-in-tree=", 22) == 0) {
      pivotLimitInTree = atoi(a + 22);
    } else if (strncmp(a, "--pivot-limit=", 14) == 0) {
      pivotLimit = atoi(a + 14);
    } else if (strncmp(a, "--away=", 7) == 0) {
      away = atof(a + 7);
    } else if (strncmp(a, "--max-cut-per-round=", 20) == 0) {
      maxCutPerRound = atoi(a + 20);
    } else if (strncmp(a, "--maximum-candidates=", 21) == 0) {
      maximumCandidates = atoi(a + 21);
    } else if (strncmp(a, "--maximum-cut-length=", 21) == 0) {
      maximumCutLength = atoi(a + 21);
    } else if (strncmp(a, "--min-violation=", 16) == 0) {
      minViolation = atof(a + 16);
    } else if (strncmp(a, "--pivot-tol=", 12) == 0) {
      pivotTol = atof(a + 12);
    } else if (strncmp(a, "--failed-pivot-limit=", 21) == 0) {
      failedPivotLimit = atoi(a + 21);
    } else if (strncmp(a, "--degenerate-pivot-limit=", 25) == 0) {
      degeneratePivotLimit = atoi(a + 25);
    } else if (strncmp(a, "--extra-cuts-limit=", 19) == 0) {
      extraCutsLimit = atoi(a + 19);
    } else if (strncmp(a, "--extra-cuts=", 13) == 0) {
      extraCuts = atoi(a + 13);
    } else if (strncmp(a, "--pivot-selection=", 18) == 0) {
      pivotSelection = atoi(a + 18);
    } else if (strncmp(a, "--sep-space=", 12) == 0) {
      sepSpace = atoi(a + 12);
    } else if (strncmp(a, "--normalization=", 16) == 0) {
      normalization = atoi(a + 16);
    } else if (strncmp(a, "--lhs-norm=", 11) == 0) {
      lhsNorm = atoi(a + 11);
    } else if (strncmp(a, "--rhs-weight-type=", 18) == 0) {
      rhsWeightType = atoi(a + 18);
    } else if (strncmp(a, "--time-limit=", 13) == 0) {
      timeLimit = atof(a + 13);
    } else if (strncmp(a, "--single-cut-time-limit=", 24) == 0) {
      singleCutTimeLimit = atof(a + 24);
    } else if (a[0] == '-') {
      fprintf(stderr, "ERROR: unknown option %s\n", a);
      usage(argv[0]);
      return 1;
    } else {
      stemArg = a;
    }
  }

  if (!stemArg) {
    fprintf(stderr, "ERROR: no fixture stem given\n");
    usage(argv[0]);
    return 1;
  }

  // Rejected rather than tolerated: without the re-solve, round 2 would separate
  // against the round-1 LP with the round-1 cuts added as rows but never priced,
  // which is not a call CBC ever makes. A single round has no such problem -- the
  // only thing lost is the bound, which the flag already disclaims.
  if (noBound && maxRounds > 1) {
    fprintf(stderr, "ERROR: --no-bound needs --rounds=1 (got %d); later rounds\n"
                    "would separate against an LP the earlier cuts never entered\n",
      maxRounds);
    return 1;
  }
  // Refused rather than passed through: :552 reads a negative maximumCutLength as
  // "already switched off" and returns immediately, so the run would report a
  // 0-cut 0-second row that looks like a spectacular speedup.
  if (maximumCutLength < 0) {
    fprintf(stderr, "ERROR: --maximum-cut-length=%d is negative, which CglLandP reads "
                    "as \"switched off\" (CglLandP.cpp:552) rather than as a tighter "
                    "cap; the run would measure an immediate return\n",
      maximumCutLength);
    return 1;
  }

  const std::string stem = fixtureStem(stemArg);

  Fixture f;
  if (!loadFixture(f, stem, quiet))
    return 1;

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();
  const int nElements0 = f.si.getNumElements();
  const int nRangeRows = countRangeRows(f.si);
  const int nCandidates = benchCandidates(f.si, away, maximumCandidates);

  int totalRowCuts = 0, totalColCuts = 0;
  size_t totalCutLen = 0;
  int maxCutLen = 0;
  double totalViol = 0.0, maxViol = 0.0;
  double totalSepTime = 0.0, totalResolveTime = 0.0;
  int totalResolveIters = 0;
  int rejSmallViol = 0, rejSmallCoef = 0, rejBigDynamic = 0, rejDense = 0,
      rejEmpty = 0;
  int rangeDisabled = 0;
  std::string rowCutsPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution is
  // reported separately: a generator whose whole gain arrives in round 1 behaves
  // differently under CBC's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A fresh generator per round, matching CBC, which constructs one per cut
    // pass. For LandP this is not merely hygiene: `numrows_` is set on the first
    // call and never reset, `extraCuts_` accumulates across calls, and the
    // range-row path negates `maximumCutLength` so a reused generator would
    // return immediately from round 2 on any range-row model.
    CglLandP landp;
    CglLandP::Parameters &p = landp.parameter();
    p.pivotLimit = pivotLimit;
    p.pivotLimitInTree = pivotLimitInTree;
    p.away = away;
    p.maxCutPerRound = maxCutPerRound;
    p.maximumCandidates = maximumCandidates;
    p.maximumCutLength = maximumCutLength;
    p.pivotTol = pivotTol;
    p.failedPivotLimit = failedPivotLimit;
    p.degeneratePivotLimit = degeneratePivotLimit;
    p.extraCutsLimit = extraCutsLimit;
    p.generateExtraCuts = (CglLandP::ExtraCutsMode)extraCuts;
    p.pivotSelection = (CglLandP::SelectionRules)pivotSelection;
    p.sepSpace = (CglLandP::SeparationSpaces)sepSpace;
    p.normalization = (CglLandP::Normalization)normalization;
    p.lhs_norm = (CglLandP::LHSnorm)lhsNorm;
    p.rhsWeightType = (CglLandP::RhsWeightType)rhsWeightType;
    p.modularize = modularize;
    p.strengthen = strengthen;
    p.perturb = perturb;
    p.countMistakenRc = countMistakenRc;
    p.exactRetry = exactRetry;
    p.exactBest = exactBest;
    p.preLengthGate = preLengthGate;
    p.useTableauRow = useTableauRow;
    p.timeLimit = timeLimit;
    p.singleCutTimeLimit = singleCutTimeLimit;
    landp.validator().setMinViolation(minViolation);
    landp.setLogLevel(logLevel);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve.
    const int nc = f.si.getNumCols();
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + nc);

    OsiCuts cs;

    // CglTreeInfo as CbcCutGenerator fills it at the root. LandP reads only
    // inTree (the pivot-limit switch) and pass (a log message).
    CglTreeInfo info;
    info.level = 0;
    info.pass = pass;
    info.formulation_rows = nRows0;
    info.inTree = inTree;
    info.options = options;

    LANDP_PROF_RESET();
    const double t0 = wallClock();
    landp.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;
    {
      // Tagged with the fixture and round so a multi-fixture serial sweep can be
      // reduced with awk without losing which row each stage belongs to.
      char tag[256];
      snprintf(tag, sizeof(tag), "%s r%d", baseName(stem).c_str(), round);
      LANDP_PROF_PRINT(tag);
    }

    // Read back rather than inferred: a negative maximumCutLength after the call
    // is the range-row self-disable having fired, and it means every later call
    // on this generator would return at once.
    if (landp.parameter().maximumCutLength < 0)
      rangeDisabled = 1;

    // Cumulative over the generator's lifetime, which for a fresh generator per
    // round is this round -- so summing across rounds is right. Undercounts
    // SmallViolation; see the file comment.
    rejSmallViol += landp.validator().numRejected(LAP::Validator::SmallViolation);
    rejSmallCoef += landp.validator().numRejected(LAP::Validator::SmallCoefficient);
    rejBigDynamic += landp.validator().numRejected(LAP::Validator::BigDynamic);
    rejDense += landp.validator().numRejected(LAP::Validator::DenseCut);
    rejEmpty += landp.validator().numRejected(LAP::Validator::EmptyCut);

    const int nRow = cs.sizeRowCuts();
    const int nCol = cs.sizeColCuts();
    if (dumpCutsFlag)
      dumpCuts(cs, baseName(stem), round);
    if (nRow == 0 && nCol == 0)
      break;

    double roundViol = 0.0;
    for (int c = 0; c < nRow; ++c) {
      const OsiRowCut &rc = cs.rowCut(c);
      const double v = rc.violated(xRound.data());
      roundViol += v;
      if (v > maxViol)
        maxViol = v;
      const int len = rc.row().getNumElements();
      totalCutLen += (size_t)len;
      // Reported alongside the mean because `maximumCutLength` caps exactly this,
      // so a sweep over it wants to see whether the cap was reached at all before
      // concluding the knob did nothing.
      if (len > maxCutLen)
        maxCutLen = len;
    }
    totalViol += roundViol;
    totalRowCuts += nRow;
    totalColCuts += nCol;

    char buf[64];
    snprintf(buf, sizeof(buf), "%s%d", round ? "+" : "", nRow);
    rowCutsPerRound += buf;

    f.si.applyCuts(cs);

    // applyCuts stays even under --no-bound: it is part of the round and costs
    // nothing next to the re-solve. Only the re-solve is skipped, so
    // isProvenOptimal still holds from the warm-started LP and the bound below
    // reads as no movement, which is what the flag documents.
    if (!noBound) {
      const double t1 = wallClock();
      f.si.resolve();
      totalResolveTime += wallClock() - t1;
      totalResolveIters += f.si.getIterationCount();
    }

    // A round whose LP does not reach optimality contributes 0 rather than a
    // bound taken from an unsolved LP.
    const double objRoundEnd
      = f.si.isProvenOptimal() ? f.si.getObjValue() : objRoundStart;
    snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "",
      f.si.getObjSense() * (objRoundEnd - objRoundStart));
    objImprovePerRound += buf;
    objRoundStart = objRoundEnd;

    ++pass;
  }

  if (rowCutsPerRound.empty()) {
    rowCutsPerRound = "0";
    objImprovePerRound = "0";
  }

  const double objEnd = f.si.isProvenOptimal() ? f.si.getObjValue() : objStart;
  const double objImprove = f.si.getObjSense() * (objEnd - objStart);
  // Relative to the starting bound's magnitude, which is what makes gains
  // comparable across instances whose objectives differ by orders of magnitude.
  const double objImproveRel
    = objImprove / (fabs(objStart) > 1.0 ? fabs(objStart) : 1.0);

  int intCols = 0;
  for (int j = 0; j < f.si.getNumCols(); ++j)
    if (!f.si.isContinuous(j))
      ++intCols;

  // Unconditional, and not silenced by --quiet: the row about to be printed
  // carries objImprove=0 and resolveTime=0, which is indistinguishable from a
  // change that destroyed the bound. Anyone reading such a row out of a log needs
  // to know why.
  if (noBound)
    fprintf(stderr, "[landp-bench] --no-bound: objEnd/objImprove/boundMoved and "
                    "resolveTime/resolveIters are NOT measured (reported as no "
                    "movement). Timing of sepTime only.\n");
  // Likewise unconditional: at pivotLimit 0 the cuts are plain MIG cuts, not
  // lift-and-project cuts, so the row is a timing control and not comparable to a
  // default row on anything but sepTime.
  if (pivotLimit == 0)
    fprintf(stderr, "[landp-bench] --pivot-limit=0: this is the generateMig() "
                    "control path -- no clones, no cut-improving pivots, DIFFERENT "
                    "CUTS. Use for timing decomposition only.\n");

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%g,%d,%d,%d,%g,%g,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g,%g,"
         "%d,%d,%d,%d,"
         "%d,%d,%d,%.10g,%.10g,%.3f,%d,"
         "%.6f,%.6f,%d,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,"
         "%d,%d,%d,%d,%d,%d,"
         "%d,%d,%d,%d,%d,%d,"
         "%d,%d,%s,%s\n",
    baseName(stem).c_str(), pivotLimit, pivotLimitInTree, away, maxCutPerRound,
    maximumCandidates, maximumCutLength, minViolation, pivotTol, failedPivotLimit,
    degeneratePivotLimit, extraCutsLimit, extraCuts, pivotSelection, sepSpace,
    normalization, lhsNorm, rhsWeightType, (int)modularize, (int)strengthen,
    (int)perturb, (int)useTableauRow, timeLimit, singleCutTimeLimit,
    pass - round, (int)inTree, options, round,
    totalRowCuts, totalColCuts, f.si.getNumRows() - nRows0, totalViol, maxViol,
    totalRowCuts ? (double)totalCutLen / totalRowCuts : 0.0, maxCutLen,
    totalSepTime, f.warmStartTime, f.warmStartIters,
    totalResolveTime, totalResolveIters,
    objStart, objEnd, objImprove, objImproveRel,
    (int)boundMoved(objStart, objImprove),
    nRows0, f.si.getNumCols(), intCols, nElements0, nRangeRows, nCandidates,
    rejSmallViol, rejSmallCoef, rejBigDynamic, rejDense, rejEmpty, rangeDisabled,
    f.restoredColTypes, (int)f.paddedRowDropped,
    rowCutsPerRound.c_str(), objImprovePerRound.c_str());

  return 0;
}
