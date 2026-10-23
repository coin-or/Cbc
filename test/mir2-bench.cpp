/**
 * Stand-alone benchmark for CglMixedIntegerRounding2.
 *
 * @file mir2-bench.cpp
 * @brief replay a CglMixedIntegerRounding2 call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcMirFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, LP solution, column types -- warm-starts the LP, and runs N
 * rounds of CglMixedIntegerRounding2. Reaching this state inside CBC costs a
 * preprocess plus a root LP, so this exists to make the iterate-measure loop on
 * the c-MIR separation itself affordable: milliseconds instead of a full solve
 * per experiment.
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
 * WHAT IS DIFFERENT ABOUT MIR2, AND WHY THIS BENCH IS NOT A COPY OF
 * twomir-bench.
 *
 * 1. `MAXAGGR_ == 1` IS CBC'S REAL CONFIGURATION, and it makes a large fraction
 *    of the generator dead. CbcSolverCutSetup.cpp:341 constructs
 *    `CglMixedIntegerRounding2 mixedGen(1, true, 1)`, and the `setMAXAGGR_` call
 *    two lines later is guarded by `if (mixedRoundStrategy != 1)` while
 *    `mixedRoundStrategy_` **defaults to 1** (CbcSolver.hpp:206, re-asserted at
 *    CbcSolver.cpp:10926). So the aggregation loop at :868 runs exactly once with
 *    `iAggregate == 0` and takes the `if` arm at :880; `selectRowToAggregate`
 *    (:1124) and `aggregateRow` (:1209) are **never called**; both SAFE_ROWS
 *    blocks inside the loop are unreachable; and the `matrixByCol` transpose
 *    built unconditionally at :195 has **no reachable reader**. Reading this
 *    generator's documented aggregation as if it ran is reading dead code.
 *    `--max-aggr=2` turns it back on -- which is not CBC's call, so it is not a
 *    measurement of CBC's cost, but it is the only way to exercise those
 *    branches, and every exactness run must include it.
 *
 * 2. THE BASIS IS THE INPUT, BUT NOT THROUGH A FACTORIZATION.
 *    `needsOptimalBasis()` is **not overridden**, so it returns false: MIR2 never
 *    factorizes and never calls getBInvARow. The `.bas` is required anyway,
 *    because `mixIntRoundPreprocess` reads `si.getRowActivity()` at :502 to
 *    resolve every range (`'R'`) row into a `'G'` or an `'L'` (:505-518). Land on
 *    a different LP point and `sense_`/`RHS_` differ, so the row types differ, so
 *    the starting-row set differs -- a different problem, not noise. This is why
 *    loadFixture refuses a fixture with no `.bas` rather than solving cold, and
 *    why `warmStartIters != 0` invalidates a comparison rather than adding to it.
 *
 * 3. IT EMITS ROW CUTS ONLY. One insertion point,
 *    `cs.insertIfNotDuplicateAndClean(cMirCut, 31, tolTest)` at :1044 with
 *    `CoinAbsFltEq tolTest(1.0e-4)` (:865). No OsiColCut anywhere, so there are
 *    no column-cut columns here and the exactness gate is simpler than Twomir's.
 *    --dump-cuts still dumps column cuts if any appear, and warns loudly: that
 *    would contradict the reading above rather than being a harmless extra.
 *
 * 4. THE ROOT CALL RUNS THE GENERATOR TWICE. `CGL_HAS_CLP` is defined
 *    (Cgl/src/CglCommon/config.h:11), so `#define MODIFY_LP 2` at :78 is live.
 *    At `!info.inTree && objSense == 1.0 && info.level >= 0` (:84-85) the branch
 *    deep-copies the whole solver (`OsiClpSolverInterface si2 = *clpSolver;`,
 *    :126), sign-flips every pure `>=` row, calls `model->setNewRowCopy(NULL)`
 *    (:150, forcing an O(nnz) row-copy rebuild), then **recurses** with
 *    `info2.level = -1-info.level` and returns (:174-176). With
 *    `setDoPreproc(1)`, `mixIntRoundPreprocess` has already run at :57-74 on the
 *    *unflipped* model and nothing between :74 and :175 reads a member it wrote,
 *    so that first preprocess pass is entirely discarded. Hence the two `--level`
 *    floors below.
 *
 *    Two details worth knowing before quoting `modifyLpFired`. `getConstClpSolver`
 *    is a plain `dynamic_cast` here -- the pointer-arithmetic variant in
 *    OsiClpSolverInterface.hpp:1661-1680 is behind `#if 0` -- so it is non-NULL
 *    for this bench's OsiClp solver, as it is for CBC's. And when `nChanged > 0`
 *    but the ClpPackedMatrix cast fails, the flip is skipped while the recursion
 *    still happens (:145-176): the recursion is outside `if (matrix)`.
 *
 * 5. `info.pass == 0` IS THE CALL WORTH MEASURING, for one reason, and it is not
 *    Twomir's. The `#if 1` block at :1001-1008 rejects a cut with
 *    `n > 0.8*numCols_` -- but only `if (info_->pass || info_->inTree)` (:1002).
 *    At the root first pass that filter is **off**, so pass 0 accepts wider cuts
 *    than any later pass. It is the analogue of Gomory's `limitAtRoot_`.
 *
 * 6. MIR2 HAS NO FIRST-PASS KILL SWITCH, the opposite of Twomir.
 *    CbcSolverCutSetup.cpp:347 passes `switches = 0 | (ALL_LAGRANGEAN *
 *    lagrangeanFlag)` with `lagrangeanFlag == 0` by default, so
 *    `setSwitchOffIfLessThan(0)` -- a condition that can never fire. A `disabled`
 *    next-run in `printGeneratorTable` therefore comes from CbcCutGenerator's own
 *    frequency tuning after the root, **not** from a zero-cuts veto on the first
 *    call. Do not attribute it to a mechanism MIR2 does not have.
 *
 * 7. NO CLOCK AND NO RANDOMNESS. Nothing in the .cpp reads CoinCpuTime,
 *    CoinWallclockTime or getMaximumSeconds, and nothing draws a random number.
 *    Replays are reproducible in a way the time-limited generators are not, and
 *    the fixture's `maxSeconds` is provenance only.
 *
 * 8. `info.formulation_rows` IS NEVER READ. The only CglTreeInfo fields MIR2
 *    touches are `inTree` (:48, :84, :215, :1002), `pass` (:48, :215, :1002),
 *    `level` (:85, :174) and `options` (:215). So there is deliberately no
 *    `--formulation-rows` flag and no `formulationRows` column here, unlike
 *    twomir-bench: offering the knob would imply it does something.
 *
 * 9. THE SETTERS THROW, and one getter is lossy -- again the opposite of Twomir,
 *    whose setAway* silently ignore out-of-range values. `setMAXAGGR_` throws
 *    CoinError unless the value is `> 0`, `-1` or `-2` (hpp:171-179);
 *    `setCRITERION_` throws outside 1..3 (hpp:191-199); `setDoPreproc` throws
 *    outside {-1, 0, 1} (cpp:1969-1978). But `getDoPreproc()` returns
 *    `doPreproc_ != 0` (cpp:1980-1983), so **-1 and 1 are indistinguishable
 *    through the getter** -- and that is exactly the distinction that decides
 *    whether preprocessing re-runs on every call. So the CSV echoes the getters
 *    for maxAggr/multiply/criterion, and for doPreproc carries the bench's own
 *    value (set inside a try/catch that hard-fails on throw) *plus*
 *    `doPreprocNonzero` from the lossy getter, so the lossiness is visible rather
 *    than papered over.
 *
 * THE CONTROL FLAGS. §7 of the process doc: price a stage before optimizing it.
 * MIR2 has no setCutTypes, so the stage gates here are call-shape flags:
 *
 *   --level=-1   the SEPARATION floor. Reproduces the inner, post-flip call: no
 *                deep copy, no sign flip, no row-copy rebuild, no recursion, one
 *                preprocess pass. This is the call that actually separates.
 *
 *   --level=0    the default, and CBC's call: preprocess, deep copy, flip,
 *                setNewRowCopy(NULL), preprocess again, then separate.
 *
 *   --max-aggr=2 turns the aggregation machinery (matrixByCol,
 *                selectRowToAggregate, aggregateRow, both SAFE_ROWS blocks) back
 *                on. Not CBC's configuration.
 *
 *   --criterion=2|3   the alternative bound-substitution rules in isLowerSubst
 *                (:1238-1251). CBC hardcodes 1 and 2/3 are never exercised.
 *
 * CAVEAT, and it is not optional when quoting a stage: the two --level runs are
 * two different calls, not two stages of one call. `level0 - level(-1)` is an
 * **attribution, not a partition** -- among other things the deep copy changes
 * which LP the separation then runs on. Attribute the accept/emit path separately
 * via rowCuts_n.
 *
 * Usage:
 *   mir2-bench <stem> [options]           (stem = dir/name.tag)
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.bas, <stem>.sol,
 * <stem>.ctype, <stem>.meta. A full path to any one of those also works; the
 * suffix is stripped.
 */

#include "CglMixedIntegerRounding2.hpp"
#include "CglTreeInfo.hpp"
#include "ClpSimplex.hpp"
#include "CoinError.hpp"
#include "CoinFinite.hpp"
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
 * Per-stage attribution, only when Cgl was built -DCGL_MIR2_PROFILE.
 *
 * Declared here rather than in CglMixedIntegerRounding2.hpp on purpose: this is a
 * local diagnostic build, not API, and putting it in the shipped header would
 * mean a header/library mismatch produces a link error for anyone who builds Cgl
 * normally. The cost of that choice is that these two prototypes must match the
 * definitions in CglMixedIntegerRounding2.cpp by hand.
 *
 * The accumulators inside the .cpp are file-static, so a profile run must be
 * serial -- which is already the rule for any run whose times get quoted. Reset
 * and print are both per round: round 0 is the call with the wide-cut filter at
 * :1002 switched off, and averaging it with the cheap later rounds would hide the
 * stage worth attacking.
 *
 * The scaffold that defines these is deliberately NOT committed. Heed its
 * lesson: instrumentation placed inside the loop being measured lands in the
 * residual, not in the stage, so time the enclosing loop as its own region and
 * require stages + instrumentation ~= loop.
 *
 * Note the MODIFY_LP==2 recursion means a --level=0 profile prints TWICE per
 * round, once per frame. Sum them, or profile at --level=-1.
 */
#ifdef CGL_MIR2_PROFILE
void cglMir2ProfileReset();
void cglMir2ProfilePrint(const char *tag);
#define MIR2_PROF_RESET() cglMir2ProfileReset()
#define MIR2_PROF_PRINT(tag) cglMir2ProfilePrint(tag)
#else
#define MIR2_PROF_RESET()
#define MIR2_PROF_PRINT(tag)
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
 * cbcMirFixtureWriteMps. For MIR2 the row is not merely cosmetic even though it
 * is redundant: it goes through determineRowType like any other row, so it can
 * add a starting row (it is all-1.0 coefficients over the empty columns, which
 * classifies by their integrality), and it would also make the `.bas` one
 * artificial short of the matrix and so unusable.
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
 * For MIR2 a lost marker does not merely shrink the problem, it reclassifies
 * rows. `mixIntRoundPreprocess` copies the whole column-type array into
 * `integerType_` at :475, and `determineRowType` (:722-810) decides ROW_MIX /
 * ROW_CONT / ROW_INT / the variable-bound types purely by counting integer
 * against continuous entries. So one lost marker can move a row out of the
 * starting set, change which rows qualify as variable bounds (`vubs_`/`vlbs_`,
 * :645-689), and change what `boundSubstitution` puts in the mixed knapsack.
 * Skipping this step measures different cuts, not a smaller problem.
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

  // getColType() caches, and mixIntRoundPreprocess copies that cached array
  // wholesale into integerType_ at :475. Without this the generator classifies
  // every row from the pre-restore view.
  si.getColType(true);
  return restored;
}

/*
 * Emit every cut, coefficient by coefficient, in a form where textual equality
 * is bit equality.
 *
 * The CSV fields are a fingerprint, not a proof: totalViol and cutNzTotal are
 * sums over all cuts, so two different cut sets can agree on every one of them.
 * When the claim being checked is "this optimization changes nothing", a
 * fingerprint match is the weaker statement. %a prints the exact IEEE double, so
 * a one-ulp difference in a single coefficient of a single cut shows up as a diff
 * rather than rounding away at the 15th digit.
 *
 * Order is as generated, deliberately unsorted, and for MIR2 that is not merely
 * convention: cuts enter through `insertIfNotDuplicateAndClean(cMirCut, 31,
 * tolTest)` with `CoinAbsFltEq(1.0e-4)` (:865, :1044), so the *sequence* decides
 * which of two near-duplicate cuts wins. A reordering is a behaviour change.
 *
 * Column cuts are dumped even though MIR2 has no OsiColCut path at all -- the
 * only insertion point in the generator is the row-cut one. If one ever appears
 * here the premise is wrong, so it is printed and warned about rather than
 * silently skipped.
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
  if (cs.sizeColCuts() > 0) {
    fprintf(stderr, "WARNING: %s: %d column cut(s) from CglMixedIntegerRounding2, "
                    "which has no OsiColCut insertion point; the bench's reading of "
                    "the generator is wrong somewhere\n",
      name.c_str(), cs.sizeColCuts());
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

/**
 * How many rows the MODIFY_LP==2 branch would sign-flip, on the LP as loaded.
 *
 * The test at :137-142 is `rowUpper[i] < 1.0e50 ? swap=false : swap=true`, i.e.
 * only rows with no finite upper bound -- pure `>=` rows. Range rows are never
 * touched, which is what keeps the `'R'` branch of preprocessing out of the
 * flip's blast radius.
 *
 * Computed live rather than read from `.meta` because it is the count that
 * decides whether *this* run takes the branch: `nChanged == 0` falls through to
 * `else { delete [] swap; }` at :177-179 and no deep copy, no flip and no
 * recursion happen.
 */
static int countGeRowsToFlip(const OsiSolverInterface &si)
{
  const double *rowUpper = si.getRowUpper();
  const int numRows = si.getNumRows();
  int n = 0;
  for (int i = 0; i < numRows; ++i)
    if (!(rowUpper[i] < 1.0e50))
      ++n;
  return n;
}

/// Everything the fixture determines, plus how long loading it took.
struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  /// Pivots the warm start needed. Zero is the expected value, and for MIR2 it is
  /// not merely a hygiene check -- see loadFixture.
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
 * either wrong looks like success while silently changing the experiment. MIR2
 * does not factorize, but the LP point is still its input: preprocessing reads
 * `si.getRowActivity()` at :502 and uses it to resolve every range row into a
 * `'G'` or an `'L'` (:505-518), so a different point gives different `sense_` and
 * `RHS_`, hence different row types, hence a different starting-row set.
 *
 * First, `readBasis` writes into ClpSimplex's `status_`, but OsiClp caches a
 * separate `CoinWarmStartBasis basis_` and `resolve()` overwrites the model from
 * it. `setWarmStart(NULL)` refreshes that cache from the model, which is what
 * makes the file's basis survive into the solve.
 *
 * Second, `resolve()` rather than `initialSolve()`: presolve discards the basis.
 * With both in place an already-optimal fixture costs 0 iterations, which is the
 * cheap self-check that the warm start worked. `setPerturbation(50)` matches the
 * generator; perturbation left on moves the vertex even from a correct basis.
 *
 * A fixture with no usable basis is refused outright rather than solved cold, and
 * that differs deliberately from probing-bench, which warns and continues. A cold
 * solve reaches *some* optimum, and for a range-row instance that is a different
 * row classification.
 *
 * One MIR2-only refusal on top of Gomory's sequence: a `.meta` reporting
 * `startingRows == 0`. The outer loop at :866 runs
 * `numRowMix_ + numRowContVB_ + numRowInt_` times, so at zero the generator
 * cannot produce anything and a timing row from it measures preprocessing alone.
 * This is MIR2's only structural precondition -- there is no instance-wide early
 * return like Twomir's free-column check, which is why the MIR2 fixture
 * population is larger than Twomir's.
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

  // Before anything else: with no starting rows the generator's outer loop runs
  // zero times, so there is no separation to measure. -1 means the key is absent
  // (a fixture written before it existed), which is not a refusal.
  const long startingRows = metaInt(meta, "startingRows", -1);
  if (startingRows == 0) {
    fprintf(stderr, "ERROR: %s: .meta reports startingRows 0, so the loop at "
                    "CglMixedIntegerRounding2.cpp:866 runs zero times and no cut can "
                    "be generated; this fixture cannot be replayed\n",
      baseName(stem).c_str());
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

  // Before the generator: integrality decides every row's type, and therefore
  // the entire starting-row set.
  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  if (!fileExists(bas)) {
    fprintf(stderr, "ERROR: %s: no .bas. MIR2 does not factorize, but preprocessing "
                    "resolves range rows from getRowActivity() "
                    "(CglMixedIntegerRounding2.cpp:502-518), so a different LP point "
                    "gives different row types and there is nothing to replay\n",
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
  // installs the stale cache over it.
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
  // means the fixture landed on a different point, which for a range-row instance
  // is a different row classification -- so the CSV carries warmStartIters and a
  // sweep should filter on it.
  if (f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured LP point "
                    "did not survive, so range rows may classify differently and these "
                    "cuts are not comparable across runs\n",
      baseName(stem).c_str(), f.warmStartIters);
  }

  f.ok = true;
  return true;
}

/**
 * Every knob, in one place, so the pre-loop probe and the per-round generator are
 * configured by the same code.
 *
 * MIR2's setters **throw** CoinError on an out-of-range value rather than
 * ignoring it (hpp:171-179, :191-199, cpp:1969-1978), so this can throw and every
 * caller must be inside a try/catch. That is the opposite hazard from Twomir,
 * whose setAway* silently keep the old value -- there, a bench echoing its own
 * variables would misreport; here it would abort. Both are handled by echoing the
 * getters, except for doPreproc, whose getter cannot represent -1 (see the file
 * comment).
 */
struct Knobs {
  // CBC's root values: CglMixedIntegerRounding2 mixedGen(1, true, 1) at
  // CbcSolverCutSetup.cpp:341 plus setDoPreproc(1) at :342. These happen to equal
  // the default constructor's gutsOfConstruct(1, true, 1, -1) except for
  // doPreproc, which CBC raises from -1 to 1.
  int maxAggr = 1;
  bool multiply = true;
  int criterion = 1;
  int doPreproc = 1;
};

static void configure(CglMixedIntegerRounding2 &g, const Knobs &k)
{
  g.setMAXAGGR_(k.maxAggr);
  g.setMULTIPLY_(k.multiply);
  g.setCRITERION_(k.criterion);
  g.setDoPreproc(k.doPreproc);
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.bas/.sol/.ctype/.meta are appended.\n"
    "Passing any one of those files also works. The .bas is REQUIRED even though\n"
    "MIR2 never factorizes: mixIntRoundPreprocess resolves every range row from\n"
    "getRowActivity() (:502-518), so a different LP point gives different row types.\n"
    "\n"
    "Call shape (defaults reproduce CBC's root call -- see the file comment):\n"
    "  --rounds=N          MIR2 rounds, LP re-solved between them (default 1; CBC\n"
    "                      calls MIR2 once per pass with a fresh generator). >1 is\n"
    "                      warned about: the wide-cut filter at :1002 turns ON as\n"
    "                      soon as info.pass is nonzero, so round 2 runs a different\n"
    "                      acceptance rule, and setDoPreproc(1) re-preprocesses\n"
    "                      against an LP the round-1 cuts have already moved.\n"
    "  --pass=N            info.pass (default 0). NOT cosmetic: pass 0 is the only\n"
    "                      pass where the `n > 0.8*numCols_` cut-width filter at\n"
    "                      :1001-1008 is switched off.\n"
    "  --in-tree           info.inTree = true. Switches the :1002 filter on, and\n"
    "                      disables the MODIFY_LP==2 branch entirely (:84).\n"
    "  --level=N           info.level (default 0). THE STAGE GATE:\n"
    "                        0  CBC's call -- preprocess, deep-copy the solver,\n"
    "                           sign-flip every >= row, setNewRowCopy(NULL),\n"
    "                           preprocess again, separate.\n"
    "                       -1  the inner post-flip call: separation only, no copy,\n"
    "                           no flip, no recursion, one preprocess pass.\n"
    "                      level0 minus level(-1) is an ATTRIBUTION, NOT A\n"
    "                      PARTITION: they are two different calls, and the deep\n"
    "                      copy changes which LP the separation then runs on.\n"
    "  --options=N         info.options bitmask (default 0, which is what\n"
    "                      CbcCutGenerator assigns at the root). Bits 4 and 8 are\n"
    "                      what :215 tests to mark cuts globally valid, so a\n"
    "                      nonzero value changes the output, not just the speed.\n"
    "                      MIR2 reads no other bit.\n"
    "\n"
    "  (There is deliberately NO --formulation-rows: MIR2 never reads\n"
    "  info.formulation_rows. The only fields it touches are inTree, pass, level\n"
    "  and options.)\n"
    "\n"
    "CBC's knobs (defaults from CbcSolverCutSetup.cpp:341-342). ALL FOUR SETTERS\n"
    "THROW CoinError on an out-of-range value rather than ignoring it, so a bad\n"
    "flag is a hard failure here, not a silent fallback:\n"
    "  --max-aggr=N        MAXAGGR_ (default 1, which is CBC's: the setMAXAGGR_ at\n"
    "                      CbcSolverCutSetup.cpp:344 is guarded by\n"
    "                      `mixedRoundStrategy != 1` and that default IS 1). At 1\n"
    "                      the aggregation loop runs once, selectRowToAggregate and\n"
    "                      aggregateRow are never called, both SAFE_ROWS blocks in\n"
    "                      the loop are unreachable and the matrixByCol transpose at\n"
    "                      :195 has no reader. Use --max-aggr=2 to exercise all of\n"
    "                      that -- required for an exactness run, but NOT a\n"
    "                      measurement of CBC's cost. 0 THROWS; -1 and -2 are\n"
    "                      accepted and mean 'up at root' (:47-56).\n"
    "  --no-multiply       MULTIPLY_ = false, which drops the inner multiplier loop\n"
    "                      bound from 2 to 1 (:838-842)\n"
    "  --criterion=N       CRITERION_ 1..3, the bound-substitution rule in\n"
    "                      isLowerSubst (:1238-1251). CBC hardcodes 1; 2 and 3 are\n"
    "                      implemented and never exercised. Outside 1..3 THROWS.\n"
    "  --do-preproc=N      doPreproc_, one of -1, 0, 1 (default 1, CBC's). Anything\n"
    "                      else THROWS. NOTE getDoPreproc() returns doPreproc_ != 0,\n"
    "                      so -1 and 1 are indistinguishable through the getter --\n"
    "                      the CSV carries this flag's own value as `doPreproc` and\n"
    "                      the getter as `doPreprocNonzero`.\n"
    "\n"
    "Output:\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --self-test         run internal consistency checks and exit\n"
    "  --dump-cuts         print every cut coefficient as an exact IEEE hex double\n"
    "                      ([cut] lines on stdout, before the CSV row). Use this,\n"
    "                      not the CSV fields, to check that a change really is\n"
    "                      output-neutral: the CSV aggregates and two different cut\n"
    "                      sets can share every aggregate. Order is as generated and\n"
    "                      is load-bearing -- insertIfNotDuplicateAndClean with\n"
    "                      CoinAbsFltEq(1.0e-4) makes the sequence decide which of\n"
    "                      two near-duplicates wins.\n"
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
 * Does `needle` appear anywhere in CglMixedIntegerRounding2.cpp?
 *
 * Used by the self-test to pin source facts the write-up depends on. Returns -1
 * when the source is not readable from the cwd, which is a skip and not a failure
 * -- the binary is installable away from the tree.
 */
static int sourceContains(const char *needle)
{
  const char *src
    = "../../Cgl/src/CglMixedIntegerRounding2/CglMixedIntegerRounding2.cpp";
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

/// Number of comma-separated fields in a CSV line.
static int countFields(const char *s)
{
  int n = 1;
  for (const char *p = s; *p; ++p)
    if (*p == ',')
      ++n;
  return n;
}

/// Number of printf conversions in a format string, ignoring "%%".
static int countConversions(const char *s)
{
  int n = 0;
  for (const char *p = s; *p; ++p) {
    if (*p != '%')
      continue;
    if (*(p + 1) == '%') {
      ++p;
      continue;
    }
    ++n;
  }
  return n;
}

/*
 * 54 fields. Gomory's and Twomir's generic set, minus everything specific to a
 * factorization or to column cuts (MIR2 has neither), plus the MIR2 axes.
 *
 * maxAggr, multiply and criterion are echoed FROM THE GETTERS, because the
 * setters throw rather than ignore and a caught-and-defaulted value must not be
 * reported as the requested one. doPreproc is the exception and is echoed from
 * the bench's own value, with the lossy getter alongside as doPreprocNonzero --
 * getDoPreproc() is `doPreproc_ != 0`, so it cannot tell -1 from 1.
 *
 * The sizing axes (startingRows ... rowNzMax) come from `.meta`, i.e. from the
 * dumper's replication of mixIntRoundPreprocess at capture time. A sweep whose
 * axis is not in the row cannot be reduced afterwards.
 *
 * A comparison script must drop the three timing columns -- sepTime,
 * warmStartTime, resolveTime -- BY NAME, never by index, leaving 51 fields whose
 * equality is checked as strings.
 */
static const char *CSV_HEADER
  = "name,maxAggr,multiply,criterion,doPreproc,doPreprocNonzero,"
    "pass,inTree,options,level,rounds,"
    "upperLimit,wideFilterActive,modifyLpFired,geRowsToFlipLive,"
    "rowCuts_n,rowsAdded,cutNzTotal,cutNzMax,avgCutLen,totalViol,maxViol,"
    "sepTime,warmStartTime,warmStartIters,resolveTime,resolveIters,"
    "objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "rows,cols,intCols,elements,"
    "startingRows,rowsMix,rowsContVB,rowsInt,rowsCont,rangeRows,geRowsToFlip,"
    "vubCount,vlbCount,cmirCalls,knapsackIntMax,cmirWork,nnzWork,rowNzMax,"
    "restoredInt,padDropped,rowCutsPerRound,objImprovePerRound";

static const char *CSV_FORMAT
  = "%s,%d,%d,%d,%d,%d,"
    "%d,%d,%d,%d,%d,"
    "%d,%d,%d,%d,"
    "%d,%d,%ld,%d,%.3f,%.10g,%.10g,"
    "%.6f,%.6f,%d,%.6f,%d,"
    "%.15g,%.15g,%.6g,%.6g,%d,"
    "%d,%d,%d,%d,"
    "%ld,%ld,%ld,%ld,%ld,%ld,%ld,"
    "%ld,%ld,%ld,%ld,%.0f,%.0f,%ld,"
    "%d,%d,%s,%s\n";

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
    { "d/x.mir.mps.gz", "d/x.mir" },
    { "d/x.mir.bas", "d/x.mir" },
    { "d/x.mir.bas.status", "d/x.mir" },
    { "d/x.mir.meta", "d/x.mir" },
    { "d/x.mir", "d/x.mir" },
  };
  for (size_t i = 0; i < sizeof(stems) / sizeof(stems[0]); ++i) {
    const std::string got = fixtureStem(stems[i].in);
    if (got != stems[i].want) {
      fprintf(stderr, "FAIL fixtureStem(%s) = %s, want %s\n",
        stems[i].in, got.c_str(), stems[i].want);
      ++failures;
    }
  }

  // The header and the row must have the same number of fields, or every
  // consumer that resolves a column BY NAME off --header reads the wrong one --
  // and reads it silently, since the types mostly happen to parse. Counting
  // conversions in the format rather than building a dummy row keeps this check
  // from needing its own duplicate argument list.
  {
    const int nHeader = countFields(CSV_HEADER);
    const int nFormat = countConversions(CSV_FORMAT);
    if (nHeader != nFormat) {
      fprintf(stderr, "FAIL CSV header has %d fields but the row format has %d "
                      "conversions; a by-name column lookup would silently read the "
                      "wrong field\n",
        nHeader, nFormat);
      ++failures;
    }
    if (strchr(CSV_HEADER, ' ')) {
      fprintf(stderr, "FAIL CSV header contains a space; the harness splits on ','\n");
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

  // Every knob the CSV reports must round-trip, because a rejected setter would
  // otherwise measure the default while the CSV reported what was asked for.
  {
    try {
      CglMixedIntegerRounding2 g;
      Knobs k;
      configure(g, k);
      if (g.getMAXAGGR_() != 1 || !g.getMULTIPLY_() || g.getCRITERION_() != 1
        || !g.getDoPreproc()) {
        fprintf(stderr, "FAIL CBC's root knobs do not round-trip: maxAggr=%d "
                        "multiply=%d criterion=%d doPreprocNonzero=%d\n",
          g.getMAXAGGR_(), (int)g.getMULTIPLY_(), g.getCRITERION_(),
          (int)g.getDoPreproc());
        ++failures;
      }
    } catch (CoinError &e) {
      fprintf(stderr, "FAIL configuring CBC's own root knobs threw: %s\n",
        e.message().c_str());
      ++failures;
    }
  }

  // The setters THROW on out-of-range rather than silently keeping the old value
  // -- the opposite of CglTwomir's setAway*. Every flag this bench offers is
  // therefore a hard-failure path, and these checks are what keep that claim
  // honest: if a setter ever starts ignoring instead, the usage text is wrong and
  // a bad flag would silently measure the default.
  {
    struct {
      int value;
      bool shouldThrow;
    } aggrs[] = { { 0, true }, { -1, false }, { -2, false }, { 7, false },
      { -3, true } };
    for (size_t i = 0; i < sizeof(aggrs) / sizeof(aggrs[0]); ++i) {
      CglMixedIntegerRounding2 g;
      bool threw = false;
      try {
        g.setMAXAGGR_(aggrs[i].value);
      } catch (CoinError &) {
        threw = true;
      }
      if (threw != aggrs[i].shouldThrow) {
        fprintf(stderr, "FAIL setMAXAGGR_(%d) %s; expected it %sto throw\n",
          aggrs[i].value, threw ? "threw" : "did not throw",
          aggrs[i].shouldThrow ? "" : "not ");
        ++failures;
      } else if (!threw && g.getMAXAGGR_() != aggrs[i].value) {
        fprintf(stderr, "FAIL setMAXAGGR_(%d) round-trips as %d\n",
          aggrs[i].value, g.getMAXAGGR_());
        ++failures;
      }
    }
  }
  {
    struct {
      int value;
      bool shouldThrow;
    } crits[] = { { 0, true }, { 1, false }, { 2, false }, { 3, false },
      { 4, true } };
    for (size_t i = 0; i < sizeof(crits) / sizeof(crits[0]); ++i) {
      CglMixedIntegerRounding2 g;
      bool threw = false;
      try {
        g.setCRITERION_(crits[i].value);
      } catch (CoinError &) {
        threw = true;
      }
      if (threw != crits[i].shouldThrow) {
        fprintf(stderr, "FAIL setCRITERION_(%d) %s; expected it %sto throw\n",
          crits[i].value, threw ? "threw" : "did not throw",
          crits[i].shouldThrow ? "" : "not ");
        ++failures;
      } else if (!threw && g.getCRITERION_() != crits[i].value) {
        fprintf(stderr, "FAIL setCRITERION_(%d) round-trips as %d\n",
          crits[i].value, g.getCRITERION_());
        ++failures;
      }
    }
  }
  {
    struct {
      int value;
      bool shouldThrow;
    } preps[] = { { -1, false }, { 0, false }, { 1, false }, { 2, true },
      { -2, true } };
    for (size_t i = 0; i < sizeof(preps) / sizeof(preps[0]); ++i) {
      CglMixedIntegerRounding2 g;
      bool threw = false;
      try {
        g.setDoPreproc(preps[i].value);
      } catch (CoinError &) {
        threw = true;
      }
      if (threw != preps[i].shouldThrow) {
        fprintf(stderr, "FAIL setDoPreproc(%d) %s; expected it %sto throw\n",
          preps[i].value, threw ? "threw" : "did not throw",
          preps[i].shouldThrow ? "" : "not ");
        ++failures;
      }
    }
  }

  // getDoPreproc() IS LOSSY: it returns doPreproc_ != 0, so -1 ("follow the
  // solver's presolve settings") and 1 ("preprocess on every call") read
  // identically. That is the one knob a bench cannot read back, which is exactly
  // why the CSV carries both `doPreproc` (this bench's value) and
  // `doPreprocNonzero` (the getter). Pinned so the two-column design does not
  // look like redundancy to a later reader.
  {
    CglMixedIntegerRounding2 g;
    g.setDoPreproc(-1);
    const bool asMinusOne = g.getDoPreproc();
    g.setDoPreproc(1);
    const bool asOne = g.getDoPreproc();
    g.setDoPreproc(0);
    const bool asZero = g.getDoPreproc();
    if (!asMinusOne || !asOne || asZero) {
      fprintf(stderr, "FAIL getDoPreproc() is documented as `doPreproc_ != 0`; got "
                      "-1->%d 1->%d 0->%d\n",
        (int)asMinusOne, (int)asOne, (int)asZero);
      ++failures;
    }
  }

  // needsOptimalBasis() is FALSE for MIR2 -- it is not overridden, and the
  // generator never factorizes. Pinned deliberately, because the .bas is required
  // here for a completely different reason (range rows resolved from
  // getRowActivity()). If this ever becomes true, the file comment's explanation
  // is no longer the whole story and needs revisiting rather than quietly aging.
  {
    CglMixedIntegerRounding2 g;
    if (g.needsOptimalBasis()) {
      fprintf(stderr, "FAIL needsOptimalBasis() is now true; this bench documents the "
                      ".bas requirement as coming from range-row resolution, not from "
                      "a factorization\n");
      ++failures;
    }
  }

  // countGeRowsToFlip must mirror the branch's own test (`rowUpper[i] < 1.0e50`),
  // including that a range row is NOT flipped -- which is the fact that bounds the
  // blast radius of anything touching that path.
  {
    OsiClpSolverInterface si;
    si.getModelPtr()->setLogLevel(0);
    const double col[1] = { 1.0 };
    const int idx[1] = { 0 };
    double lb = 0.0, ub = 1.0, obj = 1.0;
    si.addCol(0, NULL, NULL, lb, ub, obj);
    si.addRow(1, idx, col, 1.0, COIN_DBL_MAX); // pure >=, flipped
    si.addRow(1, idx, col, -COIN_DBL_MAX, 5.0); // pure <=, not flipped
    si.addRow(1, idx, col, 1.0, 5.0); // range, not flipped
    si.addRow(1, idx, col, 2.0, 2.0); // equality, not flipped
    const int got = countGeRowsToFlip(si);
    if (got != 1) {
      fprintf(stderr, "FAIL countGeRowsToFlip = %d over {G, L, R, E}, want 1\n", got);
      ++failures;
    }
  }

  // Source facts the write-up and the stage gates rest on. Silent removal of any
  // one of them invalidates the documentation, so fail loudly rather than let the
  // comments rot.
  {
    struct {
      const char *needle;
      const char *why;
    } anchors[] = {
      { "#define MODIFY_LP 2",
        "the deep-copy/sign-flip/recurse branch is what --level=0 vs --level=-1 "
        "prices; at MODIFY_LP==1 the two levels would mean something else" },
      { "indRows_[iRow] = iRow",
        "matrixByRow is documented as an identity submatrix copy of "
        "si.getMatrixByRow()" },
      { "info_->pass || info_->inTree",
        "pass 0 is documented as the only call with the 0.8*numCols_ cut-width "
        "filter switched off" },
      { "#define SAFE_ROWS",
        "the SAFE_ROWS blocks are documented as dead because MAXAGGR_ == 1, not "
        "because the macro is undefined -- if it were undefined the reader "
        "enumeration behind --max-aggr=2 would be a different list" },
    };
    for (size_t i = 0; i < sizeof(anchors) / sizeof(anchors[0]); ++i) {
      const int found = sourceContains(anchors[i].needle);
      if (found < 0) {
        printf("  (skipped source check for '%s': CglMixedIntegerRounding2.cpp not "
               "readable from cwd)\n",
          anchors[i].needle);
      } else if (!found) {
        fprintf(stderr, "FAIL '%s' is gone from CglMixedIntegerRounding2.cpp; %s\n",
          anchors[i].needle, anchors[i].why);
        ++failures;
      }
    }
  }

  printf("self-test: %d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}

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
  // One round by default: CBC constructs a fresh CglMixedIntegerRounding2 per cut
  // pass, so a single round is the call being optimized.
  int maxRounds = 1;

  Knobs k;

  int pass = 0;
  bool inTree = false;
  int options = 0;
  int level = 0;

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
    } else if (strcmp(a, "--no-multiply") == 0) {
      k.multiply = false;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--pass=", 7) == 0) {
      pass = atoi(a + 7);
    } else if (strncmp(a, "--options=", 10) == 0) {
      options = atoi(a + 10);
    } else if (strncmp(a, "--level=", 8) == 0) {
      level = atoi(a + 8);
    } else if (strncmp(a, "--max-aggr=", 11) == 0) {
      k.maxAggr = atoi(a + 11);
    } else if (strncmp(a, "--criterion=", 12) == 0) {
      k.criterion = atoi(a + 12);
    } else if (strncmp(a, "--do-preproc=", 13) == 0) {
      k.doPreproc = atoi(a + 13);
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

  const std::string stem = fixtureStem(stemArg);

  // Resolve the knobs the generator will actually run with, before the loop, so
  // the CSV can echo the getters. Unlike Twomir's silently-ignoring setters these
  // THROW, so a rejected value is a hard failure here and there is no
  // "requested vs effective" divergence to warn about -- only a bad flag to
  // report clearly instead of an uncaught CoinError and an abort.
  int effMaxAggr = 0, effCriterion = 0;
  bool effMultiply = false, effDoPreprocNonzero = false;
  try {
    CglMixedIntegerRounding2 probe;
    configure(probe, k);
    effMaxAggr = probe.getMAXAGGR_();
    effMultiply = probe.getMULTIPLY_();
    effCriterion = probe.getCRITERION_();
    effDoPreprocNonzero = probe.getDoPreproc();
  } catch (CoinError &e) {
    fprintf(stderr, "ERROR: rejected knob (--max-aggr=%d --criterion=%d "
                    "--do-preproc=%d): %s\n"
                    "       maxAggr must be > 0, -1 or -2; criterion 1..3; "
                    "doPreproc -1, 0 or 1\n",
      k.maxAggr, k.criterion, k.doPreproc, e.message().c_str());
    return 1;
  }

  Fixture f;
  if (!loadFixture(f, stem, quiet))
    return 1;

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();
  const int nElements0 = f.si.getNumElements();

  // Whether the MODIFY_LP==2 branch will actually fire, computed from the branch's
  // own conditions (:84-85, :137-145) against the LP as loaded. Without this
  // column a --level=0 measurement has no denominator: a fixture with no pure `>=`
  // row takes the `else { delete [] swap; }` at :177-179 and pays none of the
  // deep copy, the flip, the row-copy rebuild or the second preprocess pass.
  //
  // getConstClpSolver is a plain dynamic_cast (the pointer-arithmetic variant in
  // OsiClpSolverInterface.hpp is behind `#if 0`), so it is non-NULL for this
  // OsiClp solver -- reproduced here with the cast rather than the helper so the
  // bench does not depend on which variant is compiled in.
  const int geRowsLive = countGeRowsToFlip(f.si);
  const bool haveClp = dynamic_cast< const OsiClpSolverInterface * >(&f.si) != NULL;
  const bool modifyLpFired = !inTree && haveClp && f.si.getObjSense() == 1.0
    && level >= 0 && geRowsLive > 0;
  const bool wideFilterActive = (pass != 0) || inTree;
  const int upperLimit = effMultiply ? 2 : 1;

  if (!quiet && level >= 0 && !modifyLpFired) {
    fprintf(stderr, "WARNING: %s: the MODIFY_LP==2 branch will NOT fire at --level=%d "
                    "(clpSolver=%d objSense=%g inTree=%d geRowsToFlip=%d), so this run "
                    "pays no deep copy, no sign flip and one preprocess pass -- it is "
                    "not the overhead --level=0 is meant to price\n",
      baseName(stem).c_str(), level, (int)haveClp, f.si.getObjSense(), (int)inTree,
      geRowsLive);
  }
  if (!quiet && maxRounds > 1) {
    fprintf(stderr, "WARNING: %s: --rounds=%d; from round 2 on info.pass is nonzero, "
                    "which switches ON the `n > 0.8*numCols_` cut filter at "
                    "CglMixedIntegerRounding2.cpp:1002, and setDoPreproc re-runs "
                    "preprocessing against an LP the earlier cuts have moved. Later "
                    "rounds are not the call CBC makes at the root\n",
      baseName(stem).c_str(), maxRounds);
  }

  int totalRowCuts = 0;
  size_t totalCutNz = 0;
  int maxCutNz = 0;
  double totalViol = 0.0, maxViol = 0.0;
  double totalSepTime = 0.0, totalResolveTime = 0.0;
  int totalResolveIters = 0;
  std::string rowCutsPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution is
  // reported separately: a generator whose whole gain arrives in round 1 behaves
  // differently under CBC's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A fresh generator per round, matching CBC, which constructs one per cut
    // pass. It also removes any chance of state carrying between rounds and
    // confounding a sweep -- the trap that made the clique bench do the same. For
    // MIR2 the carried state would be `doneInitPre_` and every preprocessed array
    // (sense_, RHS_, rowTypes_, vubs_, vlbs_, indRows_ ...), which at
    // doPreproc_ == 1 would be recomputed anyway but at doPreproc_ != 1 would
    // silently be reused from round 1's LP.
    CglMixedIntegerRounding2 mir2;
    configure(mir2, k);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve.
    const int nc = f.si.getNumCols();
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + nc);

    OsiCuts cs;

    // The Osi-level entry, which is what CbcCutGenerator calls. Everything MIR2
    // does is inside it, including the preprocess pass and -- at level >= 0 with
    // pure `>=` rows present -- the solver deep copy and the recursive second
    // call, so all of that is inside the timed region where CBC has it.
    CglTreeInfo info;
    info.level = level;
    info.pass = pass;
    info.inTree = inTree;
    info.options = options;

    MIR2_PROF_RESET();
    const double t0 = wallClock();
    mir2.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;
    {
      // Tagged with the fixture and round so a multi-fixture serial sweep can be
      // reduced with awk without losing which row each stage belongs to.
      char tag[256];
      snprintf(tag, sizeof(tag), "%s r%d", baseName(stem).c_str(), round);
      MIR2_PROF_PRINT(tag);
    }

    const int nRow = cs.sizeRowCuts();
    if (dumpCutsFlag)
      dumpCuts(cs, baseName(stem), round);
    if (nRow == 0 && cs.sizeColCuts() == 0)
      break;

    double roundViol = 0.0;
    for (int c = 0; c < nRow; ++c) {
      const OsiRowCut &rc = cs.rowCut(c);
      const double v = rc.violated(xRound.data());
      roundViol += v;
      if (v > maxViol)
        maxViol = v;
      const int len = rc.row().getNumElements();
      totalCutNz += (size_t)len;
      // Reported alongside the mean because the `n > 0.8*numCols_` filter at
      // :1001-1008 caps exactly this -- and only when info.pass or info.inTree is
      // set, so at pass 0 a maxCutNz above that threshold is expected and is the
      // visible sign the filter was off.
      if (len > maxCutNz)
        maxCutNz = len;
    }
    totalViol += roundViol;
    totalRowCuts += nRow;

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
    fprintf(stderr, "[mir2-bench] --no-bound: objEnd/objImprove/boundMoved and "
                    "resolveTime/resolveIters are NOT measured (reported as no "
                    "movement). Timing of sepTime only.\n");

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf(CSV_FORMAT,
    baseName(stem).c_str(), effMaxAggr, (int)effMultiply, effCriterion, k.doPreproc,
    (int)effDoPreprocNonzero,
    pass - round, (int)inTree, options, level, round,
    upperLimit, (int)wideFilterActive, (int)modifyLpFired, geRowsLive,
    totalRowCuts, f.si.getNumRows() - nRows0, (long)totalCutNz, maxCutNz,
    totalRowCuts ? (double)totalCutNz / totalRowCuts : 0.0, totalViol, maxViol,
    totalSepTime, f.warmStartTime, f.warmStartIters,
    totalResolveTime, totalResolveIters,
    objStart, objEnd, objImprove, objImproveRel, (int)boundMoved(objStart, objImprove),
    nRows0, f.si.getNumCols(), intCols, nElements0,
    metaInt(meta, "startingRows", -1), metaInt(meta, "rowsMix", -1),
    metaInt(meta, "rowsContVB", -1), metaInt(meta, "rowsInt", -1),
    metaInt(meta, "rowsCont", -1), metaInt(meta, "rangeRows", -1),
    metaInt(meta, "geRowsToFlip", -1),
    metaInt(meta, "vubCount", -1), metaInt(meta, "vlbCount", -1),
    metaInt(meta, "cmirCalls", -1), metaInt(meta, "knapsackIntMax", -1),
    metaNum(meta, "cmirWork", -1.0), metaNum(meta, "nnzWork", -1.0),
    metaInt(meta, "rowNzMax", -1),
    f.restoredColTypes, (int)f.paddedRowDropped,
    rowCutsPerRound.c_str(), objImprovePerRound.c_str());

  return 0;
}
