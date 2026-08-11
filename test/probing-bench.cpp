/**
 * Stand-alone benchmark for CglProbing.
 *
 * @file probing-bench.cpp
 * @brief replay a CglProbing call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcProbingFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, LP solution, column types -- warm-starts the LP, and runs N
 * rounds of CglProbing. Reaching this state inside CBC costs a preprocess plus a
 * root LP, so this exists to make the iterate-measure loop on the probing
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
 * WHAT IS DIFFERENT ABOUT PROBING, AND WHY THIS BENCH IS NOT A COPY OF THE OTHER
 * TWO.
 *
 * 1. IT RETURNS TWO KINDS OF CUT. Row cuts (disaggregation and coefficient
 *    strengthening) go into `cs` as OsiRowCut; bound tightenings go in as
 *    OsiColCut, and *those are the ones that matter most on many instances* --
 *    they move the bound without adding a row at all. So `colCuts`,
 *    `boundsTightened` and `colCutViol` are first-class columns here, and a run
 *    that reports 0 row cuts has not necessarily done nothing. The whole-solve
 *    screen showed instances at both extremes: air05 spends 0.75s for 0 row cuts
 *    and 62 column cuts, while academictimetablesmall spends 4.93s for 6858 row
 *    cuts and 0 column cuts.
 *
 * 2. IT IS CALLED THROUGH A DIFFERENT ENTRY POINT. CbcCutGenerator.cpp:340-376
 *    special-cases CglProbing: it calls `generateCutsAndModify(si, cs, info2)`
 *    with the model's `CglTreeProbingInfo`, not `generateCuts` with a plain
 *    CglTreeInfo. The difference is not cosmetic -- generateCutsAndModify leaves
 *    colLower_/colUpper_ allocated for tightLower()/tightUpper(), and the
 *    CglTreeProbingInfo subclass makes `initializeFixing()` return 1, which turns
 *    on the saveFixingInfo path inside probe() (the `info->fixes(...)` calls that
 *    accumulate implications). A bench using generateCuts with a base CglTreeInfo
 *    would silently skip that work and measure a cheaper call than CBC's. So this
 *    bench uses generateCutsAndModify with a real CglTreeProbingInfo, and
 *    `--plain-info` exists to price the difference rather than to hide it.
 *
 * 3. `info->pass == 0` IS THE EXPENSIVE CALL, NOT MERELY THE FIRST. probe() sets
 *    `justFix = (!inTree && !pass) ? -1 : 0` and then `if (justFix < 0) maxProbe
 *    = numberThisTime_`, which *overrides* the maxProbe==123 "be a bit
 *    intelligent" subsetting that keeps only every 4th candidate. So pass 0
 *    probes every candidate column and later passes probe a quarter of them.
 *    `--pass=N` is therefore a first-class knob, and the default 0 matches both
 *    the fixture's capture point and the call worth optimizing.
 *
 * NO AUXILIARY STRUCTURE IS LOADED, AND THAT IS FAITHFUL. CglProbing does have
 * caches a fixture could have had to serialize -- rowCopy_/columnCopy_/rowLower_/
 * rowUpper_ from snapshot(), and the clique block from createCliques(). Neither
 * is populated on CBC's path: `snapshot(` appears nowhere in Cbc/src, and
 * `createCliques` appears once, on a CglKnapsackCover. So rowCopy_ is NULL,
 * gutsOfGenerateCuts rebuilds from si.getMatrixByRow() every call, numberCliques_
 * is 0 and the dispatch goes to probe() rather than probeCliques(). Rebuilding
 * from the model here is not an approximation of what CBC did -- it is exactly
 * what CBC does. `--snapshot` and `--cliques` exist to measure those paths when
 * wanted, which makes them this bench's configuration rather than a silent
 * assumption.
 *
 * THE CONTROL FLAGS. §7 of the process doc: price a stage before optimizing it.
 *
 *   --row-cuts=4     "just column cuts" -- disables disaggregation and
 *                    coefficient strengthening while leaving all the propagation
 *                    in place. time(3) - time(4) prices row-cut construction.
 *   --max-probe=0    disables probing entirely, leaving only tighten() and the
 *                    integer column cuts derived from it. This is the floor: any
 *                    fixture whose time survives --max-probe=0 is spending it in
 *                    tighten() and the minR/maxR recompute, not in probing.
 *                    Verified to be a real gate: gutsOfGenerateCuts wraps the
 *                    whole candidate-selection and probe() call in
 *                    `if (maxProbe > 0)` at CglProbing.cpp:2350.
 *
 * TWO KNOBS THAT LOOK LIKE CONTROLS AND ARE NOT. Both are exposed anyway, because
 * a sweep that silently measures a constant is the worse outcome -- but neither
 * can be used to price a stage, and the reason is recorded here so a flat result
 * is not misread as "this parameter does not matter".
 *
 *   --max-look=N     DEAD ON THIS PATH. setMaxLook/setMaxLookRoot feed
 *                    maxStack_/maxStackRoot_, and probe() reads them at
 *                    CglProbing.cpp:3616 -- then discards the value with an
 *                    unconditional `maxStack = 80;` at :3624, under a comment
 *                    reading "need a way for user to ask for more". Every one of
 *                    the ~30 later uses (stack sizing, the `nstackC < 2*maxStack`
 *                    guards in the propagation core) therefore sees 80 regardless
 *                    of what CBC or this bench asked for. There is no environment
 *                    override. Changing the stack depth requires editing the
 *                    library, so it is a code experiment, not a flag.
 *   --max-probe=N>0  does NOT bound the probing loop, despite the name. In
 *                    probe(), maxProbe is recomputed at :3975 and the
 *                    maxProbe==123 "be a bit intelligent" block only *permutes*
 *                    lookedAt_, moving the kept quarter to the front; the loop
 *                    that follows is `for (iLook = 0; iLook < numberThisTime_;
 *                    ...)` and maxProbe is never read again. The one live clamp,
 *                    `numberThisTime_ = min(numberThisTime_, maxProbe)` at :2546,
 *                    is on the snapshot branch, and its counterpart on the normal
 *                    branch (:2395) is commented out. So --max-probe is a
 *                    three-valued knob in practice: 0 disables probing, and any
 *                    positive value probes all numberThisTime_ candidates in a
 *                    slightly different order. That ordering is not nothing --
 *                    probing fixes bounds as it goes -- but it is not a budget.
 *
 * Usage:
 *   probing-bench <stem> [options]           (stem = dir/name.tag)
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.bas, <stem>.sol,
 * <stem>.ctype, <stem>.meta. A full path to any one of those also works; the
 * suffix is stripped.
 */

#include "CglProbing.hpp"
#include "CglTreeInfo.hpp"
#include "ClpSimplex.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>

/**
 * The stack depth probe() actually uses, as a literal.
 *
 * Not read from the generator, because the generator does not know it: probe()
 * overwrites its own maxStack with this constant at CglProbing.cpp:3624, after
 * reading maxStackRoot_. The CSV column exists so a --max-look sweep carries the
 * value that applied next to the value requested, making the knob's deadness
 * visible in the data rather than only in a comment. If the library's hardcode is
 * ever removed this must switch to reading the generator, so --self-test asserts
 * it against the source.
 */
static const int PROBE_HARDCODED_MAXSTACK = 80;

/// Wall-clock, not CPU: separation cost as a user would feel it.
static double wallClock()
{
  return CoinGetTimeOfDay();
}

/*
 * Per-stage attribution, only when Cgl was built -DCGL_PROBING_PROFILE.
 *
 * Declared here rather than in CglProbing.hpp on purpose: this is a local
 * diagnostic build, not API, and putting it in the shipped header would mean a
 * header/library mismatch produces a link error for anyone who builds Cgl
 * normally. The cost of that choice is that these two prototypes must match the
 * definitions in CglProbing.cpp by hand.
 *
 * The accumulators inside CglProbing.cpp are file-static, so a profile run must
 * be serial -- which is already the rule for any run whose times get quoted.
 * Reset is per round, print is per round too: a per-round breakdown is the point,
 * since round 0 (info->pass == 0) is the expensive call and averaging it with
 * the cheap later rounds would hide exactly the stage worth attacking.
 */
#ifdef CGL_PROBING_PROFILE
void cglProbingProfileReset();
void cglProbingProfilePrint(const char *tag);
#define PROBE_PROF_RESET() cglProbingProfileReset()
#define PROBE_PROF_PRINT(tag) cglProbingProfilePrint(tag)
#else
#define PROBE_PROF_RESET()
#define PROBE_PROF_PRINT(tag)
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
static long metaInt(const std::string &path, const char *key, long dflt)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return dflt;

  char line[512];
  long value = dflt;
  while (fgets(line, sizeof(line), fp)) {
    char k[256];
    double v = 0.0;
    if (sscanf(line, "%255s %lf", k, &v) == 2 && strcmp(k, key) == 0) {
      value = (long)v;
      break;
    }
  }
  fclose(fp);
  return value;
}

/**
 * Remove the pad row, if this fixture has one.
 *
 * A capture whose matrix held empty columns carries one extra final row, without
 * which writeMps would drop those columns and shift every index after them -- see
 * cbcProbeFixtureWriteMps. The row is redundant so leaving it would not change
 * the LP, but it *would* change what probing sees: gutsOfGenerateCuts counts rows
 * and filters them by length, and probe() propagates through every row a probed
 * column appears in.
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
 * For Probing that changes the experiment twice over rather than merely
 * annotating it. `intVar[i]` is the gate on whether a column enters lookedAt_ at
 * all, so lost markers shrink the candidate set and the probing loop; and it is
 * also the gate on turning a tightened bound into an OsiColCut, so lost markers
 * suppress column cuts. Skipping this step measures a smaller, easier problem.
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

  // getColType() caches, and derives Binary vs GeneralInteger from the bounds, so
  // it must be recomputed or probing would see the pre-restore view.
  si.getColType(true);
  return restored;
}

/// Everything the fixture determines, plus how long loading it took.
struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  /// Pivots the warm start needed. Zero is the expected value and the check that
  /// the captured basis actually took; see loadFixture.
  int warmStartIters = 0;
  bool paddedRowDropped = false;
  int restoredColTypes = -1;
  bool ok = false;
};

/**
 * Load the fixture and warm-start to the captured optimum.
 *
 * Two things are needed to actually *land* on the captured vertex, and getting
 * either wrong looks like success while silently changing the experiment. For
 * Probing the sensitivity is sharper than for a pure separator: gutsOfGenerateCuts
 * sorts lookedAt_ by `fabs(0.5 - frac(colsol[i]))`, so a different optimal vertex
 * reorders the probing sequence, and since probing fixes bounds as it goes, a
 * different order gives a different result, not merely a different order of the
 * same results.
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

  // Before probing: integrality gates both the candidate set and the column cuts.
  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  bool haveBasis = false;
  if (fileExists(bas)) {
    if (lp->readBasis(bas.c_str()) < 0) {
      fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
    } else {
      haveBasis = true;
      // Push the freshly-read status into OsiClp's cached basis_, or the solve
      // below installs the stale cache over it.
      f.si.setWarmStart(NULL);
    }
  } else if (!quiet) {
    fprintf(stderr, "WARNING: no basis %s; solving cold\n", bas.c_str());
  }

  lp->setPerturbation(50);

  const double t0 = wallClock();
  if (haveBasis) {
    f.si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
    f.si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
    f.si.resolve();
  } else {
    f.si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);
    f.si.initialSolve();
  }
  f.warmStartTime = wallClock() - t0;
  f.warmStartIters = f.si.getIterationCount();

  if (!f.si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: LP not optimal after warm start (%s)\n", stem.c_str());
    return false;
  }
  // A correct warm start from an optimal basis costs no pivots. Anything else
  // means the fixture landed on a different vertex, so say so rather than quietly
  // measuring a different LP -- which for Probing also means a different probing
  // order.
  if (haveBasis && f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured basis "
                    "did not survive, so this is a different vertex\n",
      baseName(stem).c_str(), f.warmStartIters);
  }

  f.ok = true;
  return true;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.bas/.sol/.ctype/.meta are appended.\n"
    "Passing any one of those files also works.\n"
    "\n"
    "Call shape (defaults reproduce CBC's root call -- see the file comment):\n"
    "  --rounds=N          probing rounds, LP re-solved between them (default 1;\n"
    "                      CBC calls Probing once per pass with a fresh generator)\n"
    "  --pass=N            info->pass (default 0). NOT cosmetic: pass 0 makes\n"
    "                      probe() override maxProbe with numberThisTime_, probing\n"
    "                      EVERY candidate; later passes keep every 4th.\n"
    "  --in-tree           info->inTree = true (selects the non-Root knob set)\n"
    "  --options=N         info->options bitmask (default 0, which is what\n"
    "                      CbcCutGenerator assigns at the root: 64/2048/16384 all\n"
    "                      clear, so probeFast and FIXED_BOTH_WAYS are unreachable)\n"
    "  --plain-info        use a base CglTreeInfo and generateCuts() instead of a\n"
    "                      CglTreeProbingInfo and generateCutsAndModify(). Prices\n"
    "                      the implication-saving path; NOT what CBC does.\n"
    "  --cutoff=F          OsiDualObjectiveLimit. With usingObjective>0 a finite\n"
    "                      cutoff appends the objective as an extra row.\n"
    "\n"
    "CBC's knobs (defaults are CBC's, from CbcSolverCutSetup.cpp:90-127, NOT\n"
    "CglProbing's constructor defaults, which differ on nearly every field):\n"
    "  --mode=N            1 = probe fractional (default), 2/3 = more\n"
    "  --max-pass=N        passes inside one call (default 1)\n"
    "  --max-probe=N       default 123. CONTROL FLAG only at 0, which disables\n"
    "                      probing and leaves tighten() plus its column cuts, so\n"
    "                      time(123) - time(0) prices probing and time(0) is the\n"
    "                      irreducible floor. A positive value does NOT cap the\n"
    "                      probing loop -- it only reorders candidates. See the\n"
    "                      file comment.\n"
    "  --max-look=N        requested stack depth (default 20). HAS NO EFFECT:\n"
    "                      probe() overwrites it with 80 unconditionally. Exposed\n"
    "                      so the CSV records requested vs actual (maxStackActual);\n"
    "                      see the file comment before sweeping it.\n"
    "  --max-elements=N    longest row considered (default 300; note\n"
    "                      gutsOfGenerateCuts resets this to nCols when the rows it\n"
    "                      would drop hold under a tenth of the nonzeros)\n"
    "  --row-cuts=N        0 none, 1 disaggregation, 2 coefficient, 3 both,\n"
    "                      4 just column cuts; -N = as +N but only fixes unless at\n"
    "                      root (default -3, CBC's). CONTROL FLAG: 4 isolates the\n"
    "                      column-cut work, so time(3) - time(4) prices row cuts.\n"
    "  --using-objective=N (default 1)\n"
    "  --max-seconds=F     generator budget (default 0 = none, which is also what\n"
    "                      CBC gives it: setMaxSeconds is applied only in\n"
    "                      CbcCutGenerator's non-Probing branch. Leave at 0 -- a\n"
    "                      wall-clock budget makes results irreproducible.)\n"
    "\n"
    "Cached-structure paths, off by default because CBC never takes them:\n"
    "  --snapshot          call snapshot() first, so rowCopy_ is populated and\n"
    "                      gutsOfGenerateCuts uses the cached matrix\n"
    "  --cliques[=min,max] call createCliques(), which routes to probeCliques()\n"
    "\n"
    "Output:\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --self-test         run internal consistency checks and exit\n"
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
    { "d/x.probe.mps.gz", "d/x.probe" },
    { "d/x.probe.bas", "d/x.probe" },
    { "d/x.probe.bas.status", "d/x.probe" },
    { "d/x.probe.meta", "d/x.probe" },
    { "d/x.probe", "d/x.probe" },
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
  // relaxation" from "current point now infeasible".
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

  // --max-look is documented as dead, and the CSV reports a literal 80 for the
  // depth that applied. Both claims rest on a line in the library, so pin them:
  // the setter must round-trip (proving the deadness is probe()'s override and not
  // a broken setter), and the hardcode must still read 80 in the source. If Cgl
  // ever honours the knob, this fails and the column stops lying rather than
  // quietly reporting a constant that is no longer true.
  {
    CglProbing p;
    p.setMaxLookRoot(37);
    if (p.getMaxLookRoot() != 37) {
      fprintf(stderr, "FAIL setMaxLookRoot/getMaxLookRoot round trip: got %d\n",
        p.getMaxLookRoot());
      ++failures;
    }
  }
  {
    // Relative to this source file, since the bench builds inside the workspace.
    const char *src = "../../Cgl/src/CglProbing/CglProbing.cpp";
    FILE *fp = fopen(src, "r");
    if (!fp) {
      printf("  (skipped maxStack hardcode check: %s not readable from cwd)\n", src);
    } else {
      char line[512];
      bool found = false;
      while (fgets(line, sizeof(line), fp)) {
        int v = 0;
        // Match the bare assignment only, not the min() clamp above it.
        if (sscanf(line, " maxStack = %d ;", &v) == 1) {
          found = true;
          if (v != PROBE_HARDCODED_MAXSTACK) {
            fprintf(stderr, "FAIL maxStack hardcode is now %d, bench reports %d\n",
              v, PROBE_HARDCODED_MAXSTACK);
            ++failures;
          }
          break;
        }
      }
      fclose(fp);
      if (!found) {
        fprintf(stderr, "FAIL no unconditional 'maxStack = N;' in %s -- the knob may "
                        "now be live, so maxStackActual must read the generator\n",
          src);
        ++failures;
      }
    }
  }

  printf("self-test: %d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}

static const char *CSV_HEADER
  = "name,mode,maxPass,maxProbe,maxLook,maxStackActual,maxElements,rowCuts,"
    "usingObjective,pass,inTree,options,plainInfo,snapshot,cliques,rounds,"
    "rowCuts_n,colCuts_n,boundsTightened,rowsAdded,totalViol,maxViol,colCutSlack,"
    "colCutViol,avgCutLen,sepTime,warmStartTime,warmStartIters,resolveTime,"
    "resolveIters,objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "numberThisTime,rows,cols,intCols,probeCandidates,probeRows,probeRowNonzeros,"
    "restoredInt,rowCutsPerRound,colCutsPerRound,objImprovePerRound";

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
  // One round by default: CBC constructs a fresh CglProbing per cut pass, so a
  // single round is the call being optimized.
  int maxRounds = 1;

  // CBC's values (CbcSolverCutSetup.cpp:90-127), not CglProbing's constructor
  // defaults. The two differ on nearly every field, and benchmarking the library
  // defaults would measure a configuration CBC never uses.
  int mode = 1;
  int maxPass = 1;
  int maxProbe = 123;
  int maxLook = 20;
  int maxElements = 300;
  int rowCuts = -3;
  int usingObjective = 1;
  double maxSeconds = 0.0;

  int pass = 0;
  bool inTree = false;
  int options = 0;
  bool plainInfo = false;
  double cutoff = 0.0;
  bool haveCutoff = false;

  bool doSnapshot = false;
  bool doCliques = false;
  int cliqueMin = 2, cliqueMax = 200;

  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (strcmp(a, "--csv-header") == 0) {
      csvHeader = true;
    } else if (strcmp(a, "--quiet") == 0) {
      quiet = true;
    } else if (strcmp(a, "--in-tree") == 0) {
      inTree = true;
    } else if (strcmp(a, "--plain-info") == 0) {
      plainInfo = true;
    } else if (strcmp(a, "--snapshot") == 0) {
      doSnapshot = true;
    } else if (strcmp(a, "--cliques") == 0) {
      doCliques = true;
    } else if (strncmp(a, "--cliques=", 10) == 0) {
      doCliques = true;
      if (sscanf(a + 10, "%d,%d", &cliqueMin, &cliqueMax) < 1) {
        fprintf(stderr, "ERROR: --cliques= wants min[,max]\n");
        return 1;
      }
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--pass=", 7) == 0) {
      pass = atoi(a + 7);
    } else if (strncmp(a, "--options=", 10) == 0) {
      options = atoi(a + 10);
    } else if (strncmp(a, "--mode=", 7) == 0) {
      mode = atoi(a + 7);
    } else if (strncmp(a, "--max-pass=", 11) == 0) {
      maxPass = atoi(a + 11);
    } else if (strncmp(a, "--max-probe=", 12) == 0) {
      maxProbe = atoi(a + 12);
    } else if (strncmp(a, "--max-look=", 11) == 0) {
      maxLook = atoi(a + 11);
    } else if (strncmp(a, "--max-elements=", 15) == 0) {
      maxElements = atoi(a + 15);
    } else if (strncmp(a, "--row-cuts=", 11) == 0) {
      rowCuts = atoi(a + 11);
    } else if (strncmp(a, "--using-objective=", 18) == 0) {
      usingObjective = atoi(a + 18);
    } else if (strncmp(a, "--max-seconds=", 14) == 0) {
      maxSeconds = atof(a + 14);
    } else if (strncmp(a, "--cutoff=", 9) == 0) {
      cutoff = atof(a + 9);
      haveCutoff = true;
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

  const std::string stem = fixtureStem(stemArg);

  Fixture f;
  if (!loadFixture(f, stem, quiet))
    return 1;

  if (haveCutoff)
    f.si.setDblParam(OsiDualObjectiveLimit, cutoff);

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();

  int totalRowCuts = 0, totalColCuts = 0, totalTightened = 0;
  size_t totalCutLen = 0;
  double totalViol = 0.0, maxViol = 0.0;
  double totalColSlack = 0.0, totalColViol = 0.0;
  double totalSepTime = 0.0, totalResolveTime = 0.0;
  int totalResolveIters = 0;
  int numberThisTime = 0;
  std::string rowCutsPerRound, colCutsPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution is
  // reported separately: a generator whose whole gain arrives in round 1 behaves
  // differently under CBC's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A fresh generator per round, matching CBC, which constructs one per cut
    // pass. It also removes any chance of state carrying between rounds and
    // confounding a sweep -- the trap that made the clique bench do the same.
    // Note totalTimesCalled_ resets with it, which matters: several blocks in
    // probe() are gated on totalTimesCalled_ == 0.
    CglProbing probing;
    probing.setMode(mode);
    probing.setUsingObjective(usingObjective);
    probing.setMaxPass(maxPass);
    probing.setMaxPassRoot(maxPass);
    probing.setMaxProbe(maxProbe);
    probing.setMaxProbeRoot(maxProbe);
    probing.setMaxLook(maxLook);
    probing.setMaxLookRoot(maxLook);
    probing.setMaxElements(maxElements);
    probing.setMaxElementsRoot(maxElements);
    probing.setRowCuts(rowCuts);
    // Left at 0 by default, which is also what CBC gives it: setMaxSeconds is
    // applied only in CbcCutGenerator's non-Probing branch, so the root Probing
    // call has no internal deadline. A budget here would make timings
    // self-referential.
    if (maxSeconds > 0.0)
      probing.setMaxSeconds(maxSeconds);

    if (doSnapshot)
      probing.snapshot(f.si, NULL, true);
    if (doCliques)
      probing.createCliques(f.si, cliqueMin, cliqueMax);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve. The bounds are kept for the
    // same reason -- a column cut is scored against the bounds it tightened.
    const int nc = f.si.getNumCols();
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + nc);
    const std::vector< double > loRound(f.si.getColLower(), f.si.getColLower() + nc);
    const std::vector< double > upRound(f.si.getColUpper(), f.si.getColUpper() + nc);

    OsiCuts cs;

    // Two entry points, deliberately: CBC uses generateCutsAndModify with a
    // CglTreeProbingInfo, whose initializeFixing() returns 1 and so switches on
    // the implication-saving path inside probe(). generateCuts with a base
    // CglTreeInfo skips it and is therefore a *cheaper call than CBC makes* --
    // which is why --plain-info exists to price it, not as the default.
    PROBE_PROF_RESET();
    const double t0 = wallClock();
    if (plainInfo) {
      CglTreeInfo info;
      info.level = 0;
      info.pass = pass;
      info.formulation_rows = nRows0;
      info.inTree = inTree;
      info.options = options;
      probing.generateCuts(f.si, cs, info);
    } else {
      CglTreeProbingInfo info(&f.si);
      info.level = 0;
      info.pass = pass;
      info.formulation_rows = nRows0;
      info.inTree = inTree;
      info.options = options;
      info.hasParent = 0;
      info.parentSolver = NULL;
      probing.generateCutsAndModify(f.si, cs, &info);
    }
    totalSepTime += wallClock() - t0;
    {
      // Tagged with the fixture and round so a multi-fixture serial sweep can be
      // reduced with awk without losing which row each stage belongs to.
      char tag[256];
      snprintf(tag, sizeof(tag), "%s r%d", baseName(stem).c_str(), round);
      PROBE_PROF_PRINT(tag);
    }

    numberThisTime = probing.numberThisTime();

    const int nRow = cs.sizeRowCuts();
    const int nCol = cs.sizeColCuts();
    if (nRow == 0 && nCol == 0)
      break;

    double roundViol = 0.0;
    for (int c = 0; c < nRow; ++c) {
      const OsiRowCut &rc = cs.rowCut(c);
      const double v = rc.violated(xRound.data());
      roundViol += v;
      if (v > maxViol)
        maxViol = v;
      totalCutLen += (size_t)rc.row().getNumElements();
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

    const double t1 = wallClock();
    f.si.resolve();
    totalResolveTime += wallClock() - t1;
    totalResolveIters += f.si.getIterationCount();

    // A round whose LP does not reach optimality contributes 0 rather than a
    // bound taken from an unsolved LP. Note an infeasible LP is a real outcome
    // here and not a failure: probing can prove infeasibility outright.
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

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
         "%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.3f,%.6f,%.6f,%d,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,%d,%d,%d,%d,%ld,%ld,%ld,%d,%s,%s,%s\n",
    baseName(stem).c_str(), mode, maxPass, maxProbe, maxLook,
    PROBE_HARDCODED_MAXSTACK,
    maxElements, rowCuts, usingObjective, pass - round, (int)inTree, options,
    (int)plainInfo, (int)doSnapshot, (int)doCliques, round,
    totalRowCuts, totalColCuts, totalTightened, f.si.getNumRows() - nRows0,
    totalViol, maxViol, totalColSlack, totalColViol,
    totalRowCuts ? (double)totalCutLen / totalRowCuts : 0.0,
    totalSepTime, f.warmStartTime, f.warmStartIters,
    totalResolveTime, totalResolveIters, objStart, objEnd, objImprove,
    objImproveRel, (int)boundMoved(objStart, objImprove), numberThisTime,
    nRows0, f.si.getNumCols(), intCols,
    metaInt(meta, "probeCandidates", -1), metaInt(meta, "probeRows", -1),
    metaInt(meta, "probeRowNonzeros", -1),
    f.restoredColTypes,
    rowCutsPerRound.c_str(), colCutsPerRound.c_str(), objImprovePerRound.c_str());

  return 0;
}
