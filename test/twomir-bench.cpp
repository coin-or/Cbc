/**
 * Stand-alone benchmark for CglTwomir.
 *
 * @file twomir-bench.cpp
 * @brief replay a CglTwomir call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcTwomirFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, LP solution, column types -- warm-starts the LP, and runs N
 * rounds of CglTwomir. Reaching this state inside CBC costs a preprocess plus a
 * root LP, so this exists to make the iterate-measure loop on the two-step MIR
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
 * WHAT IS DIFFERENT ABOUT TWOMIR, AND WHY THIS BENCH IS NOT A COPY OF
 * gomory-bench.
 *
 * 1. THE BASIS IS THE ALGORITHM'S INPUT, as for Gomory, but through a different
 *    mechanism worth knowing when reading a diff. needsOptimalBasis() returns
 *    true (CglTwomir.cpp:2241-2246) and DGG_generateTabRowCuts builds its own
 *    stack-local CoinFactorization from the warm start (:1564-1569), does one
 *    updateColumnTranspose per fractional basic integer (:1044), then forms the
 *    tableau row by dotting that BTRAN result against the raw column-major matrix
 *    (:1056-1067). It never calls getBInvARow. Two different optimal bases of the
 *    same degenerate vertex give entirely different cuts, both valid, so
 *    warmStartIters != 0 does not merely add noise -- it invalidates the
 *    comparison. The CSV carries the column so a sweep can be filtered on it.
 *
 * 2. A FREE COLUMN DISQUALIFIES THE WHOLE INSTANCE. CglTwomir.cpp:95-118 counts
 *    columns with lb < -1e20 && ub > 1e20 and, if any exist, **returns with no
 *    cuts at all** before doing any work. This is the single reason the Gomory or
 *    Probing fixture population cannot be reused here, and it is why loadFixture
 *    refuses such a fixture outright: replaying one measures the free-column
 *    scan and nothing else.
 *
 * 3. IT RETURNS TWO KINDS OF CUT, and the reason differs from Gomory's. A cut
 *    whose support is a single column is turned into an OsiColCut (:455-460), or,
 *    when the implied bounds cross, into a deliberately infeasible 1<=0 row cut
 *    (:443-448). So `colCuts_n`, `boundsTightened` and `colCutViol` are
 *    first-class columns, a run reporting 0 row cuts has not necessarily done
 *    nothing, and any exactness gate on a change to this generator must compare
 *    *both* kinds. Twomir's singleton branch is exactly where a bound-only
 *    behaviour change would hide.
 *
 * 4. FOUR INDEPENDENTLY GATEABLE STAGES, which is strictly more control than
 *    Gomory offered. See "THE CONTROL FLAGS" below; the short version is that
 *    setCutTypes gives *two nested floors* rather than one, and the outer one
 *    does not even build a factorization.
 *
 * 5. `info.pass == 0` IS THE CALL WORTH MEASURING, for three reasons, none of
 *    them Gomory's:
 *      - `max_elements` becomes `useSolver->getNumCols()` at :300-302 -- not the
 *        250 CBC set, not the 50000 default. The cut-length cap is at its
 *        loosest exactly here.
 *      - the tableau stage is alive only while `info.level < 1 && info.pass < 6`
 *        (:314). Pass 6 and later is a *structurally different* call that skips
 *        DGG_generateTabRowCuts entirely, not a cheaper version of the same one.
 *      - CBC gives TwoMirCuts `switches = 1`, i.e. setSwitchOffIfLessThan(1)
 *        (CbcSolverCutSetup.cpp:383). Returning zero row cuts on the first root
 *        pass disables the generator for the whole solve. So pass 0 is the only
 *        call that always happens, and the population it defines is exactly the
 *        population that kill switch judges.
 *
 * 6. NO CLOCK IS READ ANYWHERE IN CglTwomir -- no CoinCpuTime, no
 *    CoinWallclockTime, no maxSeconds. Replays are reproducible in a way the
 *    time-limited generators are not, and the fixture's `maxSeconds` is
 *    provenance only.
 *
 * 7. THE RNG IS SEEDED IN THE CONSTRUCTOR, NOT PASSED IN. `randomNumberGenerator_`
 *    is built with 987654321 (:569) under COIN_OWN_RANDOM_32, and
 *    `info.randomNumberGenerator` is ignored. Only the formulation stage draws
 *    from it (:317-320). A fresh generator per round therefore restarts the LCG,
 *    which makes each round reproducible but means round 2 is not a continuation
 *    of the stream CBC would have had. Single-round is the only configuration
 *    comparable to CBC.
 *
 * 8. WITH twomirType_ == 0 THE CLP PROLOGUE IS DEAD, so this bench calls
 *    generateCuts directly and still measures CBC's path. The
 *    CGL_HAS_CLP_TWOMIR block at :97-281 is gated on `clpSolver`, which is NULL
 *    unless `originalSolver_` is non-NULL, which only passInOriginalSolver makes
 *    true -- and CBC does that only under -latwomir (CbcSolverCutSetup.cpp:402).
 *    So on the default root path `useSolver` stays `&si`, no objective is
 *    perturbed, and the post-pass that erases locally-useless cuts (:497-530) is
 *    skipped. `--orig-solver` exposes the path deliberately; see its caveats.
 *
 * THE CONTROL FLAGS. §7 of the process doc: price a stage before optimizing it.
 * setCutTypes(mir, twomir, tab, form) gives two *nested* floors:
 *
 *   --no-tab --no-form   the PER-CALL FLOOR. DGG_getData (one getWarmStart() plus
 *                        the cast, the scratch block :748-763, the O(ncol) column
 *                        loop, the O(elements) row dot for slack activity and the
 *                        O(elements) integer-slack test), an empty emit loop, and
 *                        DGG_freeData. **No CoinFactorization is built at all** --
 *                        it is a local of DGG_generateTabRowCuts (:1564-1569).
 *                        O(elements + ncol + nrow), and no change to the inner
 *                        loops can beat it.
 *
 *   --no-mir --no-2mir   the BASE-CONSTRUCTION FLOOR. These do not disable base
 *                        construction: they set t_max = t_min-1 and
 *                        q_max = q_min-1 (:311-312) so the t- and q-loop bodies
 *                        never run, while the factorization, every
 *                        DGG_getTableauConstraint, DGG_transformConstraint and
 *                        DGG_nicefyConstraint still do. This isolates exactly the
 *                        cost of building tableau rows.
 *
 *   --no-form            prices the tableau stage.
 *   --no-tab             prices the formulation stage.
 *
 * CAVEAT, and it is not optional when quoting a stage: `floor + tab + form` only
 * accounts for the total up to the *shared emit loop* (:332-488), whose cost
 * scales with `cut_list.n`, which these flags also change. Attribute the emit
 * loop separately via rowCuts_n/colCuts_n; never fold it into a stage.
 *
 * A KNOB THAT LOOKS LIKE A CONTROL AND IS WEAKER THAN IT LOOKS.
 * `--away-at-root=0.5` zeroes tableau candidates (the test at :1578 is
 * `frac < threshold || frac > 1-threshold`) but does **not** stop the formulation
 * stage, which needs no fractional basic integer. And values above 0.5 are
 * *silently ignored* by setAwayAtRoot (:2258-2263) rather than clamped, so the
 * CSV echoes the getters, never the requested value. Use --no-tab for the
 * tableau stage.
 *
 * TWO HARD WIDTH CAPS ARE PAID AFTER THE WORK, which matters when reading a cut
 * count against a time. `base->nz > 500` (:1592) discards a tableau row that has
 * already been fully built, and `c->nz > 500` in DGG_isCutDesirable (:2174) runs
 * after DGG_substituteSlacks has already densified the cut. Time spent on cuts
 * that are then thrown away is real work with no output.
 *
 * ONE UNCONDITIONAL printf. `printf("2mir_test: why does constraint not exist ?")`
 * at :1588 is not behind `talk` or any ifdef. Any harness comparing stdout
 * between two builds must filter it -- and must report a change in how many times
 * it fires, because that IS a behaviour change.
 *
 * Usage:
 *   twomir-bench <stem> [options]           (stem = dir/name.tag)
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.bas, <stem>.sol,
 * <stem>.ctype, <stem>.meta. A full path to any one of those also works; the
 * suffix is stripped.
 */

#include "CglTwomir.hpp"
#include "CglTreeInfo.hpp"
#include "ClpSimplex.hpp"
#include "CoinTime.hpp"
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

/*
 * Per-stage attribution, only when Cgl was built -DCGL_TWOMIR_PROFILE.
 *
 * Declared here rather than in CglTwomir.hpp on purpose: this is a local
 * diagnostic build, not API, and putting it in the shipped header would mean a
 * header/library mismatch produces a link error for anyone who builds Cgl
 * normally. The cost of that choice is that these two prototypes must match the
 * definitions in CglTwomir.cpp by hand.
 *
 * The accumulators inside CglTwomir.cpp are file-static, so a profile run must be
 * serial -- which is already the rule for any run whose times get quoted. Reset
 * and print are both per round: round 0 (info.pass == 0) is the expensive call
 * with max_elements at getNumCols() and the tableau stage alive, and averaging it
 * with the cheap later rounds would hide exactly the stage worth attacking.
 *
 * The scaffold that defines these is deliberately NOT committed. Heed its
 * lesson: instrumentation placed inside the loop being measured lands in the
 * residual, not in the stage, so time the enclosing loop as its own region and
 * require stages + instrumentation ~= loop.
 */
#ifdef CGL_TWOMIR_PROFILE
void cglTwomirProfileReset();
void cglTwomirProfilePrint(const char *tag);
#define TWO_PROF_RESET() cglTwomirProfileReset()
#define TWO_PROF_PRINT(tag) cglTwomirProfilePrint(tag)
#else
#define TWO_PROF_RESET()
#define TWO_PROF_PRINT(tag)
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
 * cbcTwomirFixtureWriteMps. For Twomir the row is not merely cosmetic even though
 * it is redundant: it is an extra basis position, so leaving it in changes the
 * factorization every tableau row is derived from, and it would also make the
 * `.bas` one artificial short of the matrix and so unusable.
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
 * For Twomir a lost marker changes the experiment in four separate places, all in
 * CglTwomir.cpp. `si.getColType()` is read at :330 to build the `intVar` the emit
 * loop's singleton branch tests at :428, so a lost marker turns an OsiColCut into
 * nothing at all. DGG_getData reads it at :761 to set DGG_setIsInteger, which
 * gates tableau candidate selection at :1576 -- and, worse, uses it at :804-806 to
 * *tighten* the stored bounds with ceil/floor, so data->lb/ub differ and every
 * cut derived from them differs. It gates the integer-slack classification, which
 * decides whether a slack may carry an integral coefficient. Skipping this step
 * does not merely measure a smaller problem, it measures different cuts.
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

  // getColType() caches, and CglTwomir reads it through si.getColType() at :330
  // AND inside DGG_getData at :761, so it must be recomputed or the generator
  // sees the pre-restore view in both places.
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
 * order and later rounds depend on it. Note CglTwomir only ever sorts the *slack
 * tail* of a tableau row (CoinSort_2 at :1099 starts at nz1), never the cut list.
 *
 * Column cuts are dumped for a Twomir-specific reason, not for symmetry: the
 * singleton branch at :455-460 is the only place a bound change can appear, so a
 * comparison of row cuts alone would miss it entirely.
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
  /// Pivots the warm start needed. Zero is the expected value, and for Twomir it
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
 * Twomir, as for Gomory, the basis is the algorithm's input: landing on a
 * different optimal basis of the same vertex produces a different, equally valid
 * set of cuts. A before/after comparison across such a load is meaningless, not
 * merely noisy.
 *
 * First, `readBasis` writes into ClpSimplex's `status_`, but OsiClp caches a
 * separate `CoinWarmStartBasis basis_` and `resolve()` overwrites the model from
 * it. `setWarmStart(NULL)` refreshes that cache from the model, which is what
 * makes the file's basis survive into the solve. It is also what makes
 * `si.getWarmStart()` -- which is where DGG_getData reads the basis at :710 --
 * return the captured one.
 *
 * Second, `resolve()` rather than `initialSolve()`: presolve discards the basis.
 * With both in place an already-optimal fixture costs 0 iterations, which is the
 * cheap self-check that the warm start worked. `setPerturbation(50)` matches the
 * generator; perturbation left on moves the vertex even from a correct basis.
 *
 * A fixture with no usable basis is refused outright rather than solved cold, and
 * that differs deliberately from probing-bench, which warns and continues. A cold
 * solve reaches *some* optimal basis, which for Probing is a slightly different
 * candidate order and for Twomir is a different problem.
 *
 * Two Twomir-only refusals on top of Gomory's sequence:
 *
 *   - a `.meta` reporting freeColumns > 0 is refused. generateCuts returns at
 *     :118 having produced nothing, so replaying it times the free-column scan
 *     and reports zero cuts -- a row that looks like a catastrophic regression
 *     and is a property of the instance.
 *   - negativeSlackRows > 0 is warned about loudly. The dumper computes it from
 *     the same activity dot DGG_getData does, and at a genuine optimum no row is
 *     more than 1e-5 outside its bound, so a nonzero value means the *loader*
 *     reconstructed a different LP than the one captured.
 */
static bool loadFixture(Fixture &f, const std::string &stem, bool quiet)
{
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";
  const std::string meta = stem + ".meta";

  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no problem file for stem %s\n", stem.c_str());
    return false;
  }

  // Before anything else: an instance with a free column produces no cuts at all,
  // so there is no separation to measure and a timing row from it would be a
  // measurement of CglTwomir.cpp:95-118 and nothing more.
  const long freeColumns = metaInt(meta, "freeColumns", 0);
  if (freeColumns > 0) {
    fprintf(stderr, "ERROR: %s: .meta reports %ld free column(s). CglTwomir returns "
                    "at CglTwomir.cpp:118 without generating anything, so this "
                    "fixture cannot be replayed\n",
      baseName(stem).c_str(), freeColumns);
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

  // Before the generator: integrality gates candidate selection, tightens the
  // bounds DGG_getData stores, and decides whether a singleton becomes a column
  // cut at all.
  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  if (!fileExists(bas)) {
    fprintf(stderr, "ERROR: %s: no .bas. A two-step MIR cut is derived from a row of "
                    "the simplex tableau at a specific basis, so there is nothing to "
                    "replay without it\n",
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
  // installs the stale cache over it -- and getWarmStart() would hand DGG_getData
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
  // same one, which for Twomir is a different algorithm input -- so the CSV
  // carries warmStartIters and a sweep should filter on it.
  if (f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured basis "
                    "did not survive, so these are cuts from a DIFFERENT basis and "
                    "not comparable across runs\n",
      baseName(stem).c_str(), f.warmStartIters);
  }

  // The dumper should never have written a fixture with a violated row, so this
  // says the reconstruction differs from the capture, not that the capture was
  // odd.
  const long negativeSlackRows = metaInt(meta, "negativeSlackRows", 0);
  if (negativeSlackRows > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: .meta reports %ld row(s) with activity more than "
                    "1e-5 outside their bound at capture; the reconstructed LP is not "
                    "the captured one\n",
      baseName(stem).c_str(), negativeSlackRows);
  }

  f.ok = true;
  return true;
}

/**
 * Every knob, in one place, so the pre-loop probe and the per-round generator are
 * configured by the same code.
 *
 * This exists because CglTwomir's setters are not uniformly faithful:
 * setAway/setAwayAtRoot silently ignore a value outside (0, 0.5] rather than
 * clamping it (:2246-2263). A bench that echoed its own variables would print
 * what was requested while the generator ran something else -- the exact shape of
 * a false result. So the CSV echoes the getters, and `configure` is called on a
 * throwaway generator first so those getters can be read even before round 1.
 */
struct Knobs {
  // CBC's root values (CbcSolverCutSetup.cpp:359-418), NOT CglTwomir's
  // constructor defaults, which differ: the constructor uses away 0.0005 both at
  // the root and in the tree, and max_elements 50000.
  int maxElements = 250;
  int maxElementsRoot = 50000;
  double awayAtRoot = 0.005;
  double away = 0.01;
  int tMin = 1, tMax = 1, qMin = 1, qMax = 1, aMax = 2;
  bool doMir = true, do2mir = true, doTab = true, doForm = true;
  int twomirType = 0;
};

static void configure(CglTwomir &g, const Knobs &k)
{
  g.setMaxElements(k.maxElements);
  g.setMaxElementsRoot(k.maxElementsRoot);
  g.setAwayAtRoot(k.awayAtRoot);
  g.setAway(k.away);
  g.setMirScale(k.tMin, k.tMax);
  g.setTwomirScale(k.qMin, k.qMax);
  g.setAMax(k.aMax);
  g.setCutTypes(k.doMir, k.do2mir, k.doTab, k.doForm);
  g.setTwomirType(k.twomirType);
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.bas/.sol/.ctype/.meta are appended.\n"
    "Passing any one of those files also works. The .bas is REQUIRED: a two-step MIR\n"
    "cut is derived from a row of the simplex tableau at a specific basis.\n"
    "\n"
    "Call shape (defaults reproduce CBC's root call -- see the file comment):\n"
    "  --rounds=N          Twomir rounds, LP re-solved between them (default 1;\n"
    "                      CBC calls Twomir once per pass with a fresh generator).\n"
    "                      >1 is warned about: the tableau stage dies at pass 6\n"
    "                      (CglTwomir.cpp:314) and a fresh generator restarts the\n"
    "                      RNG the formulation stage draws from, so later rounds\n"
    "                      are not the calls CBC makes.\n"
    "  --pass=N            info.pass (default 0). NOT cosmetic: pass 0 makes\n"
    "                      max_elements getNumCols() (:300-302) rather than the 250\n"
    "                      CBC set, AND the tableau stage runs only while pass < 6.\n"
    "  --in-tree           info.inTree = true, which switches max_elements from\n"
    "                      maxElementsRoot to maxElements and the fractionality\n"
    "                      threshold from awayAtRoot to away (:296-297)\n"
    "  --level=N           info.level (default 0). level >= 1 kills the tableau\n"
    "                      stage outright (:314).\n"
    "  --options=N         info.options bitmask (default 0, which is what\n"
    "                      CbcCutGenerator assigns at the root). Bit 32 makes\n"
    "                      max_elements getNumCols() even at pass > 0 (:300); bits\n"
    "                      4 and 8 make cuts globally valid (:490-497). A nonzero\n"
    "                      value changes the output, not just the speed.\n"
    "  --formulation-rows=N   info.formulation_rows, which caps the formulation\n"
    "                      stage's outer loop (:1616). Default: the fixture's row\n"
    "                      count, which is what CBC passes at the root.\n"
    "\n"
    "CBC's knobs (defaults are CBC's, from CbcSolverCutSetup.cpp:359-418, NOT\n"
    "CglTwomir's constructor defaults, which differ 20x on away):\n"
    "  --max-elements=N    cut length cap in the tree (default 250)\n"
    "  --max-elements-root=N  cap at the root (default 50000; note CBC never calls\n"
    "                      setMaxElementsRoot at all, so 50000 is the constructor\n"
    "                      value, and at pass 0 it is overridden by getNumCols())\n"
    "  --away-at-root=F    root fractionality threshold (default 0.005, which is\n"
    "                      CBC's MORE_CUTS value -- and MORE_CUTS is #defined inside\n"
    "                      the *Gomory* block at CbcSolverCutSetup.cpp:139 and never\n"
    "                      undef'd, so this default depends on another generator's\n"
    "                      setup). Values outside (0,0.5] are SILENTLY IGNORED, so\n"
    "                      read the awayAtRoot CSV column, not your command line.\n"
    "                      Weaker as a control than --no-tab: it zeroes tableau\n"
    "                      candidates but leaves the formulation stage running.\n"
    "  --away=F            in-tree threshold (default 0.01; same silent refusal)\n"
    "  --t-min=N --t-max=N   MIR scalings (default 1,1)\n"
    "  --q-min=N --q-max=N   two-step MIR scalings (default 1,1)\n"
    "  --a-max=N           maximum bhat/alpha (default 2)\n"
    "  --twomir-type=N     0 normal (default). 1/2 (and the 10+/20+ variants) need\n"
    "                      an original solver to do anything, so use --orig-solver.\n"
    "  --orig-solver       passInOriginalSolver(clone of the fixture solver), which\n"
    "                      is what makes --twomir-type live (it sets twomirType_=1\n"
    "                      itself if it is 0, CglTwomir.cpp:656-657). CBC's default\n"
    "                      root does NOT do this -- only -latwomir does\n"
    "                      (CbcSolverCutSetup.cpp:402). Refused with --rounds>1: the\n"
    "                      CLP path relaxes free-column bounds in originalSolver_\n"
    "                      and never restores them (:148-160), and skips the\n"
    "                      objective restore when the simplex fails (:274-280), so\n"
    "                      round 2 would run against a mutated solver.\n"
    "\n"
    "Stage gates -- setCutTypes, giving two NESTED floors (see the file comment):\n"
    "  --no-tab --no-form  the per-call floor. DGG_getData + an empty emit loop +\n"
    "                      DGG_freeData. No CoinFactorization is built at all.\n"
    "  --no-mir --no-2mir  the base-construction floor. Sets t_max=t_min-1 and\n"
    "                      q_max=q_min-1 (:311-312), so every tableau row is still\n"
    "                      built and transformed but no cut body runs.\n"
    "  --no-form           prices the tableau stage.\n"
    "  --no-tab            prices the formulation stage.\n"
    "                      CAVEAT: floor+tab+form accounts for the total only up to\n"
    "                      the shared emit loop (:332-488), whose cost scales with\n"
    "                      the number of cuts, which these flags also change.\n"
    "                      Attribute it separately via rowCuts_n/colCuts_n.\n"
    "\n"
    "Output:\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --self-test         run internal consistency checks and exit\n"
    "  --dump-cuts         print every cut coefficient as an exact IEEE hex\n"
    "                      double ([cut] lines on stdout, before the CSV row), row\n"
    "                      cuts AND column cuts. Use this, not the CSV fields, to\n"
    "                      check that a change really is output-neutral: the CSV\n"
    "                      aggregates and two different cut sets can share every\n"
    "                      aggregate. NOTE CglTwomir.cpp:1588 prints an\n"
    "                      unconditional '2mir_test: why does constraint not exist'\n"
    "                      line, so a diff of two stdout streams must filter it --\n"
    "                      and must report a change in how often it fires.\n"
    "  --no-bound          skip the post-separation LP re-solve. TIMING AID ONLY:\n"
    "                      it makes objEnd/objImprove/boundMoved meaningless (they\n"
    "                      report no movement), so never judge a change with it --\n"
    "                      the bound is the metric that decides. Refused with\n"
    "                      --rounds>1, where round 2 would otherwise separate\n"
    "                      against an LP the round-1 cuts never entered.\n"
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
 * How much a column cut actually tightens, and whether it cuts off the current
 * point.
 *
 * There is no OsiColCut::violated(), and the row-cut notion does not transfer: a
 * column cut is violated when the incumbent LP value falls outside the new bound.
 * Both numbers are reported because they answer different questions -- `slack`
 * says the relaxation got smaller, `viol` says the current vertex is now
 * infeasible and so the next resolve must move.
 */
static void scoreColCut(const OsiColCut &cc, const double *x,
  const double *colLower, const double *colUpper,
  int &nTightened, double &slack, double &viol)
{
  const CoinPackedVector &lbs = cc.lbs();
  const CoinPackedVector &ubs = cc.ubs();

  const int nl = lbs.getNumElements();
  const int *li = lbs.getIndices();
  const double *lv = lbs.getElements();
  for (int k = 0; k < nl; ++k) {
    const int j = li[k];
    if (lv[k] > colLower[j] + 1.0e-12) {
      ++nTightened;
      slack += lv[k] - colLower[j];
      if (lv[k] > x[j] + 1.0e-7)
        viol += lv[k] - x[j];
    }
  }

  const int nu = ubs.getNumElements();
  const int *ui = ubs.getIndices();
  const double *uv = ubs.getElements();
  for (int k = 0; k < nu; ++k) {
    const int j = ui[k];
    if (uv[k] < colUpper[j] - 1.0e-12) {
      ++nTightened;
      slack += colUpper[j] - uv[k];
      if (uv[k] < x[j] - 1.0e-7)
        viol += x[j] - uv[k];
    }
  }
}

/**
 * Does `needle` appear anywhere in CglTwomir.cpp?
 *
 * Used by the self-test to pin two source facts the write-up depends on. Returns
 * -1 when the source is not readable from the cwd, which is a skip and not a
 * failure -- the binary is installable away from the tree.
 */
static int sourceContains(const char *needle)
{
  const char *src = "../../Cgl/src/CglTwomir/CglTwomir.cpp";
  FILE *fp = fopen(src, "r");
  if (!fp)
    return -1;
  char line[1024];
  int found = 0;
  while (fgets(line, sizeof(line), fp)) {
    if (strstr(line, needle)) {
      found = 1;
      break;
    }
  }
  fclose(fp);
  return found;
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
  // directory name that happens to end in one of them.
  struct {
    const char *in;
    const char *want;
  } stems[] = {
    { "d/x.twomir.mps.gz", "d/x.twomir" },
    { "d/x.twomir.bas", "d/x.twomir" },
    { "d/x.twomir.bas.status", "d/x.twomir" },
    { "d/x.twomir.meta", "d/x.twomir" },
    { "d/x.twomir", "d/x.twomir" },
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

  // scoreColCut must count only genuine tightenings, and must separate "smaller
  // relaxation" from "current point now infeasible". This matters more here than
  // in probing-bench: CglTwomir reaches OsiColCut only through the singleton
  // branch at CglTwomir.cpp:455-460, so a miscount there would read as
  // "generated nothing" on precisely the fixtures where a bound-only change
  // would hide.
  {
    const double x[2] = { 0.5, 0.5 };
    const double lo[2] = { 0.0, 0.0 };
    const double up[2] = { 1.0, 1.0 };
    OsiColCut cc;
    CoinPackedVector ubs;
    ubs.insert(0, 0.0); // real tightening, and cuts off x[0]=0.5
    ubs.insert(1, 1.0); // no tightening at all
    cc.setUbs(ubs);
    int n = 0;
    double slack = 0.0, viol = 0.0;
    scoreColCut(cc, x, lo, up, n, slack, viol);
    if (n != 1 || fabs(slack - 1.0) > 1.0e-12 || fabs(viol - 0.5) > 1.0e-12) {
      fprintf(stderr, "FAIL scoreColCut: n=%d slack=%g viol=%g, want 1 1 0.5\n",
        n, slack, viol);
      ++failures;
    }
  }

  // Every knob the CSV reports must round-trip, because a rejected setter would
  // otherwise measure the default while the CSV reported what was asked for.
  {
    CglTwomir g;
    Knobs k;
    configure(g, k);
    if (g.getMaxElements() != 250 || g.getMaxElementsRoot() != 50000
      || fabs(g.getAwayAtRoot() - 0.005) > 1.0e-15
      || fabs(g.getAway() - 0.01) > 1.0e-15
      || g.getTmin() != 1 || g.getTmax() != 1 || g.getQmin() != 1
      || g.getQmax() != 1 || g.getAmax() != 2 || g.twomirType() != 0) {
      fprintf(stderr, "FAIL CBC's root knobs do not round-trip: maxElements=%d "
                      "maxElementsRoot=%d awayAtRoot=%g away=%g t=%d,%d q=%d,%d "
                      "aMax=%d type=%d\n",
        g.getMaxElements(), g.getMaxElementsRoot(), g.getAwayAtRoot(), g.getAway(),
        g.getTmin(), g.getTmax(), g.getQmin(), g.getQmax(), g.getAmax(),
        g.twomirType());
      ++failures;
    }
  }

  // setAway/setAwayAtRoot REFUSE out-of-range values rather than clamping them
  // (CglTwomir.cpp:2246-2263). The CSV therefore has to echo the getters, and
  // this is the check that keeps that reasoning honest: if the setters ever start
  // clamping, echoing the getter is still right but the docs are wrong.
  {
    CglTwomir g;
    g.setAwayAtRoot(0.005);
    g.setAway(0.01);
    g.setAwayAtRoot(0.6);
    if (fabs(g.getAwayAtRoot() - 0.005) > 1.0e-15) {
      fprintf(stderr, "FAIL setAwayAtRoot(0.6) changed the value to %g; it is "
                      "documented as silently ignored, not clamped\n",
        g.getAwayAtRoot());
      ++failures;
    }
    g.setAwayAtRoot(0.0);
    if (fabs(g.getAwayAtRoot() - 0.005) > 1.0e-15) {
      fprintf(stderr, "FAIL setAwayAtRoot(0.0) changed the value to %g\n",
        g.getAwayAtRoot());
      ++failures;
    }
    g.setAway(0.6);
    if (fabs(g.getAway() - 0.01) > 1.0e-15) {
      fprintf(stderr, "FAIL setAway(0.6) changed the value to %g\n", g.getAway());
      ++failures;
    }
    // 0.5 is in range and must be accepted, or --away-at-root=0.5 would silently
    // leave 0.005 in place and "the weak control" would be a full run.
    g.setAwayAtRoot(0.5);
    if (fabs(g.getAwayAtRoot() - 0.5) > 1.0e-15) {
      fprintf(stderr, "FAIL setAwayAtRoot(0.5) rejected; got %g\n", g.getAwayAtRoot());
      ++failures;
    }
  }

  // The four stage gates are the control flags this bench's stage accounting is
  // built on, so all 16 combinations must round-trip through the getters. A gate
  // that silently did not stick would make a floor read as a full run.
  {
    for (int mask = 0; mask < 16; ++mask) {
      const bool mir = (mask & 1) != 0;
      const bool twomir = (mask & 2) != 0;
      const bool tab = (mask & 4) != 0;
      const bool form = (mask & 8) != 0;
      CglTwomir g;
      g.setCutTypes(mir, twomir, tab, form);
      if ((g.getIfMir() != 0) != mir || (g.getIfTwomir() != 0) != twomir
        || (g.getIfTableau() != 0) != tab || (g.getIfFormulation() != 0) != form) {
        fprintf(stderr, "FAIL setCutTypes(%d,%d,%d,%d) round-trips as (%d,%d,%d,%d)\n",
          (int)mir, (int)twomir, (int)tab, (int)form, g.getIfMir(), g.getIfTwomir(),
          g.getIfTableau(), g.getIfFormulation());
        ++failures;
      }
    }
  }

  // setFormulationRows is DEAD: form_nrows_ is written by the constructor, the
  // copy constructor, the assignment operator and this setter, and read nowhere.
  // The formulation stage's row cap comes from info.formulation_rows (:317-319),
  // which is why this bench exposes --formulation-rows and not the setter. Pinned
  // so nobody later mistakes one for the other -- the check is that the setter
  // does not disturb anything observable.
  {
    CglTwomir g;
    g.setMaxElements(250);
    g.setFormulationRows(7);
    if (g.getMaxElements() != 250 || g.getIfFormulation() == 0) {
      fprintf(stderr, "FAIL setFormulationRows(7) perturbed observable state; the "
                      "bench assumes form_nrows_ is write-only\n");
      ++failures;
    }
  }

  // needsOptimalBasis() is the premise of this whole bench -- it is why the .bas
  // is required rather than optional, and why a nonzero warmStartIters
  // invalidates a comparison. Pin it: if it ever returns false, that reasoning
  // needs revisiting rather than silently persisting in the comments.
  {
    CglTwomir g;
    if (!g.needsOptimalBasis()) {
      fprintf(stderr, "FAIL needsOptimalBasis() is now false; this bench's insistence "
                      "on the captured basis needs rethinking\n");
      ++failures;
    }
  }

  // Two source facts the write-up rests on. `#undef DGG_DEBUG_DGG` at :31 is what
  // makes every #if DGG_DEBUG_DGG block dead despite the header defining it 1, so
  // a note calling DGG_isConstraintViolated dead depends on it. The hard 500-width
  // cap is paid AFTER the tableau row has been fully built, which is the premise
  // of the early-exit idea and of any statement about work that produces no
  // output. Silent removal of either invalidates the documentation, so fail loudly
  // rather than let the comments rot.
  {
    struct {
      const char *needle;
      const char *why;
    } anchors[] = {
      { "#undef DGG_DEBUG_DGG",
        "every #if DGG_DEBUG_DGG block is documented as dead" },
      { "if (base->nz > 500) continue;",
        "the 500-width cap is documented as paid after the row is built" },
    };
    for (size_t i = 0; i < sizeof(anchors) / sizeof(anchors[0]); ++i) {
      const int found = sourceContains(anchors[i].needle);
      if (found < 0) {
        printf("  (skipped source check for '%s': CglTwomir.cpp not readable "
               "from cwd)\n",
          anchors[i].needle);
      } else if (!found) {
        fprintf(stderr, "FAIL '%s' is gone from CglTwomir.cpp; %s\n",
          anchors[i].needle, anchors[i].why);
        ++failures;
      }
    }
  }

  printf("self-test: %d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}

/*
 * 61 fields. Gomory's 50, minus the three factorization knobs that have no
 * CglTwomir analogue (altFactorization, condMultiplier,
 * largestFactorMultiplier), plus fourteen: the five scaling knobs (t/q/aMax) and
 * the four stage gates, without which a sweep's axis is not in the row;
 * maxElements *and* maxElementsRoot, so a reader can see that pass 0 overrode
 * both with getNumCols(); level and formulationRows, which are call-shape inputs
 * CBC supplies and which turn stages off; and formBaseRows/formIntRows, without
 * which a formulation-stage regression has no denominator.
 *
 * awayAtRoot and away are echoed FROM THE GETTERS, because the setters silently
 * ignore out-of-range values.
 *
 * A comparison script must drop the three timing columns -- sepTime,
 * warmStartTime, resolveTime -- BY NAME, never by index, leaving 58 fields whose
 * equality is checked as strings.
 */
static const char *CSV_HEADER
  = "name,maxElements,maxElementsRoot,awayAtRoot,away,tMin,tMax,qMin,qMax,aMax,"
    "doMir,do2mir,doTab,doForm,twomirType,origSolver,pass,inTree,options,level,"
    "formulationRows,rounds,"
    "rowCuts_n,colCuts_n,boundsTightened,rowsAdded,totalViol,maxViol,colCutSlack,"
    "colCutViol,avgCutLen,maxCutLen,sepTime,warmStartTime,warmStartIters,"
    "resolveTime,resolveIters,objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "rows,cols,intCols,elements,twomirCandidates,basicIntegers,nonbasicStructural,"
    "nonbasicStructuralNz,tabRowWork,outerLoopWork,activeSlackRows,"
    "integerSlackRows,formBaseRows,formIntRows,restoredInt,padDropped,"
    "rowCutsPerRound,colCutsPerRound,objImprovePerRound";

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
  // One round by default: CBC constructs a fresh CglTwomir per cut pass, so a
  // single round is the call being optimized.
  int maxRounds = 1;

  Knobs k;
  bool origSolver = false;

  int pass = 0;
  bool inTree = false;
  int options = 0;
  int level = 0;
  // -1 means "the fixture's row count", which is what CBC passes at the root
  // (numberRowsAtContinuous_). Resolved after the fixture loads.
  int formulationRows = -1;

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
    } else if (strcmp(a, "--orig-solver") == 0) {
      origSolver = true;
    } else if (strcmp(a, "--no-mir") == 0) {
      k.doMir = false;
    } else if (strcmp(a, "--no-2mir") == 0) {
      k.do2mir = false;
    } else if (strcmp(a, "--no-tab") == 0) {
      k.doTab = false;
    } else if (strcmp(a, "--no-form") == 0) {
      k.doForm = false;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--pass=", 7) == 0) {
      pass = atoi(a + 7);
    } else if (strncmp(a, "--options=", 10) == 0) {
      options = atoi(a + 10);
    } else if (strncmp(a, "--level=", 8) == 0) {
      level = atoi(a + 8);
    } else if (strncmp(a, "--formulation-rows=", 19) == 0) {
      formulationRows = atoi(a + 19);
    } else if (strncmp(a, "--max-elements-root=", 20) == 0) {
      k.maxElementsRoot = atoi(a + 20);
    } else if (strncmp(a, "--max-elements=", 15) == 0) {
      k.maxElements = atoi(a + 15);
    } else if (strncmp(a, "--away-at-root=", 15) == 0) {
      k.awayAtRoot = atof(a + 15);
    } else if (strncmp(a, "--away=", 7) == 0) {
      k.away = atof(a + 7);
    } else if (strncmp(a, "--t-min=", 8) == 0) {
      k.tMin = atoi(a + 8);
    } else if (strncmp(a, "--t-max=", 8) == 0) {
      k.tMax = atoi(a + 8);
    } else if (strncmp(a, "--q-min=", 8) == 0) {
      k.qMin = atoi(a + 8);
    } else if (strncmp(a, "--q-max=", 8) == 0) {
      k.qMax = atoi(a + 8);
    } else if (strncmp(a, "--a-max=", 8) == 0) {
      k.aMax = atoi(a + 8);
    } else if (strncmp(a, "--twomir-type=", 14) == 0) {
      k.twomirType = atoi(a + 14);
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
  // Also rejected: the CLP path mutates originalSolver_ and does not always undo
  // it. It relaxes free-column bounds to +/-1e10 and never restores them
  // (CglTwomir.cpp:148-160), and when the resolve fails it drops the saved
  // objective without putting it back (:274-280). Round 2 would then run against
  // a different problem than round 1, silently.
  if (origSolver && maxRounds > 1) {
    fprintf(stderr, "ERROR: --orig-solver needs --rounds=1 (got %d); the CLP path\n"
                    "mutates the cloned solver's bounds and objective without always\n"
                    "restoring them, so later rounds run against a different LP\n",
      maxRounds);
    return 1;
  }

  const std::string stem = fixtureStem(stemArg);

  Fixture f;
  if (!loadFixture(f, stem, quiet))
    return 1;

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();
  const int nElements0 = f.si.getNumElements();
  if (formulationRows < 0)
    formulationRows = nRows0;

  // Resolve the knobs the generator will actually run with, before the loop, so
  // the CSV can echo the getters even though the setters may have refused. Runs
  // on a throwaway object; the per-round generators are configured identically.
  double effAwayAtRoot = 0.0, effAway = 0.0;
  {
    CglTwomir probe;
    configure(probe, k);
    effAwayAtRoot = probe.getAwayAtRoot();
    effAway = probe.getAway();
    if (!quiet && fabs(effAwayAtRoot - k.awayAtRoot) > 1.0e-15) {
      fprintf(stderr, "WARNING: --away-at-root=%g was REFUSED (outside (0,0.5]); the "
                      "generator runs %g\n",
        k.awayAtRoot, effAwayAtRoot);
    }
    if (!quiet && fabs(effAway - k.away) > 1.0e-15) {
      fprintf(stderr, "WARNING: --away=%g was REFUSED (outside (0,0.5]); the generator "
                      "runs %g\n",
        k.away, effAway);
    }
  }

  int totalRowCuts = 0, totalColCuts = 0, totalTightened = 0;
  size_t totalCutLen = 0;
  int maxCutLen = 0;
  double totalViol = 0.0, maxViol = 0.0;
  double totalColSlack = 0.0, totalColViol = 0.0;
  double totalSepTime = 0.0, totalResolveTime = 0.0;
  int totalResolveIters = 0;
  std::string rowCutsPerRound, colCutsPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution is
  // reported separately: a generator whose whole gain arrives in round 1 behaves
  // differently under CBC's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A round at pass >= 6 is a structurally different call, not a cheaper one:
    // DGG_generateTabRowCuts is skipped entirely (CglTwomir.cpp:314). Warn once,
    // at the boundary, rather than let a multi-round sweep silently average a
    // tableau-stage call with formulation-only calls.
    if (pass == 6 && k.doTab && !quiet) {
      fprintf(stderr, "WARNING: %s: round %d runs at pass 6, where the tableau stage "
                      "is switched off by CglTwomir.cpp:314. Later rounds measure the "
                      "formulation stage only\n",
        baseName(stem).c_str(), round);
    }

    // A fresh generator per round, matching CBC, which constructs one per cut
    // pass. It also removes any chance of state carrying between rounds and
    // confounding a sweep -- the trap that made the clique bench do the same. For
    // Twomir the carried state would be originalSolver_ and the RNG; note the
    // RNG is reseeded to 987654321 by this construction, so each round's
    // formulation stage draws the same stream, which is reproducible but is not
    // the stream CBC's later passes would have had.
    CglTwomir twomir;
    configure(twomir, k);
    // Clones internally, so this is a per-round cost and deliberately inside the
    // loop rather than hoisted -- CBC would pay it per pass too, under -latwomir.
    // Note it sets twomirType_ = 1 itself when the type is still 0
    // (CglTwomir.cpp:656-657), so --orig-solver alone is enough to reach the CLP
    // path.
    if (origSolver)
      twomir.passInOriginalSolver(&f.si);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve. The bounds are kept for the
    // same reason -- a column cut is scored against the bounds it tightened.
    const int nc = f.si.getNumCols();
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + nc);
    const std::vector< double > loRound(f.si.getColLower(), f.si.getColLower() + nc);
    const std::vector< double > upRound(f.si.getColUpper(), f.si.getColUpper() + nc);

    OsiCuts cs;

    // The Osi-level entry, which is what CbcCutGenerator calls. With
    // twomirType_ == 0 and no original solver the CLP prologue is skipped
    // (clpSolver is NULL), so this is CBC's path, and the getColType() call at
    // :330 and the getWarmStart() inside DGG_getData stay inside the timed region
    // where CBC has them.
    CglTreeInfo info;
    info.level = level;
    info.pass = pass;
    info.formulation_rows = formulationRows;
    info.inTree = inTree;
    info.options = options;

    TWO_PROF_RESET();
    const double t0 = wallClock();
    twomir.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;
    {
      // Tagged with the fixture and round so a multi-fixture serial sweep can be
      // reduced with awk without losing which row each stage belongs to.
      char tag[256];
      snprintf(tag, sizeof(tag), "%s r%d", baseName(stem).c_str(), round);
      TWO_PROF_PRINT(tag);
    }

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
      // Reported alongside the mean because max_elements caps exactly this, so a
      // sweep over --max-elements wants to see whether the cap was reached at all
      // before concluding the knob did nothing. At pass 0 the effective cap is
      // getNumCols(), so a maxCutLen well below the flag is expected.
      if (len > maxCutLen)
        maxCutLen = len;
    }
    for (int c = 0; c < nCol; ++c) {
      scoreColCut(cs.colCut(c), xRound.data(), loRound.data(), upRound.data(),
        totalTightened, totalColSlack, totalColViol);
    }
    totalViol += roundViol;
    totalRowCuts += nRow;
    totalColCuts += nCol;

    char buf[64];
    snprintf(buf, sizeof(buf), "%s%d", round ? "+" : "", nRow);
    rowCutsPerRound += buf;
    snprintf(buf, sizeof(buf), "%s%d", round ? "+" : "", nCol);
    colCutsPerRound += buf;

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
    // bound taken from an unsolved LP. Note infeasibility is a real outcome here
    // and not a failure: the singleton branch emits a 1<=0 row cut precisely to
    // say the node is dead (CglTwomir.cpp:443-448).
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
    colCutsPerRound = "0";
    objImprovePerRound = "0";
  }

  const double objEnd = f.si.isProvenOptimal() ? f.si.getObjValue() : objStart;
  const double objImprove = f.si.getObjSense() * (objEnd - objStart);
  // Relative to the starting bound's magnitude, which is what makes gains
  // comparable across instances whose objectives differ by orders of magnitude.
  const double objImproveRel
    = objImprove / (fabs(objStart) > 1.0 ? fabs(objStart) : 1.0);

  const std::string meta = stem + ".meta";
  int intCols = 0;
  for (int j = 0; j < f.si.getNumCols(); ++j)
    if (!f.si.isContinuous(j))
      ++intCols;

  // Unconditional, and not silenced by --quiet: the row about to be printed
  // carries objImprove=0 and resolveTime=0, which is indistinguishable from a
  // change that destroyed the bound. Anyone reading such a row out of a log needs
  // to know why.
  if (noBound)
    fprintf(stderr, "[twomir-bench] --no-bound: objEnd/objImprove/boundMoved and "
                    "resolveTime/resolveIters are NOT measured (reported as no "
                    "movement). Timing of sepTime only.\n");

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%g,%g,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
         "%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.3f,%d,%.6f,%.6f,%d,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%.0f,%.0f,%ld,%ld,"
         "%ld,%ld,%d,%d,%s,%s,%s\n",
    baseName(stem).c_str(), k.maxElements, k.maxElementsRoot, effAwayAtRoot,
    effAway, k.tMin, k.tMax, k.qMin, k.qMax, k.aMax,
    (int)k.doMir, (int)k.do2mir, (int)k.doTab, (int)k.doForm, k.twomirType,
    (int)origSolver, pass - round, (int)inTree, options, level, formulationRows,
    round,
    totalRowCuts, totalColCuts, totalTightened, f.si.getNumRows() - nRows0,
    totalViol, maxViol, totalColSlack, totalColViol,
    totalRowCuts ? (double)totalCutLen / totalRowCuts : 0.0, maxCutLen,
    totalSepTime, f.warmStartTime, f.warmStartIters,
    totalResolveTime, totalResolveIters, objStart, objEnd, objImprove,
    objImproveRel, (int)boundMoved(objStart, objImprove),
    nRows0, f.si.getNumCols(), intCols, nElements0,
    metaInt(meta, "twomirCandidates", -1), metaInt(meta, "basicIntegers", -1),
    metaInt(meta, "nonbasicStructural", -1),
    metaInt(meta, "nonbasicStructuralNz", -1),
    metaNum(meta, "tabRowWork", -1.0), metaNum(meta, "outerLoopWork", -1.0),
    metaInt(meta, "activeSlackRows", -1), metaInt(meta, "integerSlackRows", -1),
    metaInt(meta, "formBaseRows", -1), metaInt(meta, "formIntRows", -1),
    f.restoredColTypes, (int)f.paddedRowDropped,
    rowCutsPerRound.c_str(), colCutsPerRound.c_str(), objImprovePerRound.c_str());

  return 0;
}
