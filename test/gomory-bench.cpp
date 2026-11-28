/**
 * Stand-alone benchmark for CglGomory.
 *
 * @file gomory-bench.cpp
 * @brief replay a CglGomory call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcGomoryFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, LP solution, column types -- warm-starts the LP, and runs N
 * rounds of CglGomory. Reaching this state inside CBC costs a preprocess plus a
 * root LP, so this exists to make the iterate-measure loop on the Gomory
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
 * WHAT IS DIFFERENT ABOUT GOMORY, AND WHY THIS BENCH IS NOT A COPY OF THE OTHER
 * THREE.
 *
 * 1. THE BASIS IS THE ALGORITHM'S INPUT, NOT PROVENANCE. needsOptimalBasis()
 *    returns true, and a Gomory cut *is* a row of the simplex tableau at the
 *    basis it is handed: CglGomory factorizes the basis it reads out of the
 *    solver's warm start (CglGomory.cpp:645-685) and BTRANs a unit vector
 *    through it. Two different optimal bases of the same degenerate vertex give
 *    entirely different cuts, both valid. For Probing a lost basis meant a
 *    reordered candidate list; here it means a different algorithm input. So
 *    warmStartIters != 0 is not a warning to note and move past -- it
 *    invalidates the comparison, and the CSV carries the column so a sweep can
 *    be filtered on it.
 *
 * 2. IT RETURNS TWO KINDS OF CUT, for a reason unrelated to Probing's. A Gomory
 *    cut over a single column is not a row cut but a bound: CglGomory.cpp:
 *    1731-1793 converts a singleton into an OsiColCut (or, when the implied
 *    bounds cross, into a deliberately infeasible 1<=0 row cut signalling that
 *    the node is dead). So `colCuts_n`, `boundsTightened` and `colCutViol` are
 *    first-class columns, a run reporting 0 row cuts has not necessarily done
 *    nothing, and any exactness gate on a change to this generator must compare
 *    *both* kinds -- comparing only row cuts would miss the bound changes
 *    entirely.
 *
 * 3. `info.pass == 0` IS THE EXPENSIVE CALL FOR TWO INDEPENDENT REASONS, and one
 *    of them is not a budget. CglGomory.cpp:842-875:
 *      - `limit = limitAtRoot_` (CBC sets 1000, or 2000 above 5000 columns)
 *        rather than `limit_` (50). That is the cap on cut length, so pass 0
 *        admits cuts 20-40x longer, and length is where the work is.
 *      - pass 0 also sets `tolerance1 = 1.0`, which at :1246 makes
 *        `useTolerance` the plain 1e-7 rather than a value scaled by the cut's
 *        own magnitude. The accuracy test at :1252 is therefore an *absolute*
 *        one at pass 0 and a relative one later: a different acceptance rule,
 *        not merely a different budget. A change validated only at pass 50 has
 *        not been validated on CBC's root call.
 *    Hence `--pass=N` is a first-class knob and the default 0 matches both the
 *    fixture's capture point and the call worth optimizing.
 *
 * 4. THERE IS NO CACHED AUXILIARY STRUCTURE TO SERIALIZE, and unlike Probing that
 *    is not a coincidence to be checked per call site -- CglGomory has nowhere to
 *    put one. Its only mutable state is `originalSolver_` (a clone, NULL unless
 *    passInOriginalSolver is called, which CBC's root setup does not do) and
 *    `numberTimesStalled_`, which is dead: both of its readers are commented out
 *    (:855, :864). Everything else is rebuilt per call from the matrix. So the
 *    fixture payload is exactly the five generic files, and rebuilding here is
 *    not an approximation of what CBC did -- it is what CBC does.
 *
 * 5. WITH gomoryType_ == 0 THE Osi-LEVEL ENTRY IS A PASS-THROUGH, so this bench
 *    calls it directly and still measures CBC's path. The CGL_HAS_CLP_GOMORY
 *    prologue at :101-274 is gated on `clpSolver`, which is NULL when
 *    `originalSolver_` is NULL, so `useSolver` stays `&si`, no objective is
 *    perturbed, and the post-pass that erases locally-useless cuts (:291-319) is
 *    skipped. `--gomory-type=N` exposes the other modes, but note they need
 *    passInOriginalSolver to do anything at all, so the flag alone changes
 *    nothing -- it is here to price `--orig-solver`, not on its own.
 *
 * THE CONTROL FLAG. §7 of the process doc: price a stage before optimizing it.
 *
 *   --away-at-root=0.5   the floor. `away` gates candidate selection at :888-901
 *                        (`reducedValue < 1-away && reducedValue > away`), so
 *                        0.5 admits nothing but an exactly-half-integral column
 *                        and drives nCandidates to ~0. What remains timed is the
 *                        irreducible per-call cost: the basis factorization
 *                        (:661-677, which is O(basis) and cannot be avoided), the
 *                        dense rowActivity build (:728-738) and the rowType
 *                        classification (:747-795). So
 *                        time(0.005) - time(0.5) prices the entire candidate loop
 *                        including the dense tableau-row dot product at
 *                        :1037-1115 that is the primary suspect, and time(0.5) is
 *                        the floor no change to the loop can beat.
 *
 * A KNOB THAT LOOKS LIKE A CONTROL AND IS NOT. `--limit-at-root=0` does *not*
 * disable anything: :852 turns 0 into `numberColumns` and :876 then adds
 * `numberRows` on top, so 0 is the *loosest* setting, not the tightest. Use a
 * small positive value to restrict cut length. Both are exposed, because a sweep
 * that silently measures the opposite of what it names is the worse outcome.
 *
 * ONE ACCEPTANCE PATH IS UNREACHABLE AT CBC's ROOT, which matters when reading a
 * cut count. Under MORE_GOMORY_CUTS==1 (the shipped value) a cut that fails the
 * accuracy test but passes a second one goes to `secondaryCuts` (:1797-1806), and
 * that container is only flushed into `cs` inside the `doSorted` block. `doSorted`
 * is `(info.options & 256) != 0` (:816) and CBC's root passes options == 0. So at
 * the default those cuts are computed and thrown away: real work that produces no
 * output. Anything that makes them cheaper is a legitimate saving, and anything
 * that makes them *appear* is a behaviour change, not a speedup.
 *
 * Usage:
 *   gomory-bench <stem> [options]           (stem = dir/name.tag)
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.bas, <stem>.sol,
 * <stem>.ctype, <stem>.meta. A full path to any one of those also works; the
 * suffix is stripped.
 */

#include "CglGomory.hpp"
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

/// Wall-clock, not CPU: separation cost as a user would feel it.
static double wallClock()
{
  return CoinGetTimeOfDay();
}

/*
 * Per-stage attribution, only when Cgl was built -DCGL_GOMORY_PROFILE.
 *
 * Declared here rather than in CglGomory.hpp on purpose: this is a local
 * diagnostic build, not API, and putting it in the shipped header would mean a
 * header/library mismatch produces a link error for anyone who builds Cgl
 * normally. The cost of that choice is that these two prototypes must match the
 * definitions in CglGomory.cpp by hand.
 *
 * The accumulators inside CglGomory.cpp are file-static, so a profile run must be
 * serial -- which is already the rule for any run whose times get quoted. Reset
 * and print are both per round: round 0 (info.pass == 0) is the expensive call
 * with the loose length limit and the absolute accuracy test, and averaging it
 * with the cheap later rounds would hide exactly the stage worth attacking.
 */
#ifdef CGL_GOMORY_PROFILE
void cglGomoryProfileReset();
void cglGomoryProfilePrint(const char *tag);
#define GOM_PROF_RESET() cglGomoryProfileReset()
#define GOM_PROF_PRINT(tag) cglGomoryProfilePrint(tag)
#else
#define GOM_PROF_RESET()
#define GOM_PROF_PRINT(tag)
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
 * cbcGomoryFixtureWriteMps. For Gomory the row is not merely cosmetic even though
 * it is redundant: it is an extra basis position, so leaving it in changes the
 * factorization the cuts are derived from, and it would also make the `.bas` one
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
 * For Gomory a lost marker changes the experiment in three separate places, all
 * in CglGomory.cpp: `intVar[]` gates candidate selection at :892; it selects the
 * mixed-integer coefficient formula at :1078-1104 rather than the continuous one,
 * *and* a column taking the continuous branch there increments
 * `numberNonInteger`, which changes the rhs relaxation ladder at :1606-1638 --
 * so the cut's right-hand side moves even for cuts that still get generated; and
 * it gates the `rowType |= 4` integer-slack test at :780, which decides whether a
 * slack variable may contribute an integral coefficient. Skipping this step does
 * not merely measure a smaller problem, it measures different cuts.
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

  // getColType() caches, and CglGomory reads it through si.getColType() at :70 to
  // build intVar[], so it must be recomputed or the generator sees the
  // pre-restore view.
  si.getColType(true);
  return restored;
}

/// Everything the fixture determines, plus how long loading it took.
/*
 * Emit every cut, coefficient by coefficient, in a form where textual equality
 * is bit equality.
 *
 * The 48 CSV fields are a fingerprint, not a proof: totalViol and avgCutLen are
 * sums over all cuts, so two different cut sets can agree on every one of them.
 * When the claim being checked is "this optimization changes nothing", a
 * fingerprint match is the weaker statement. %a prints the exact IEEE double, so
 * a one-ulp difference in a single coefficient of a single cut shows up as a
 * diff rather than rounding away at the 15th digit.
 *
 * Order is as generated, deliberately unsorted: the order cuts land in OsiCuts
 * is itself part of what must not change, since applyCuts consumes them in that
 * order and later rounds depend on it.
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

struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  /// Pivots the warm start needed. Zero is the expected value, and for Gomory it
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
 * Gomory the stakes are higher than for any other generator here: the basis is
 * the algorithm's input, so landing on a different optimal basis of the same
 * vertex produces a different, equally valid set of cuts. A before/after
 * comparison across such a load is meaningless, not merely noisy.
 *
 * First, `readBasis` writes into ClpSimplex's `status_`, but OsiClp caches a
 * separate `CoinWarmStartBasis basis_` and `resolve()` overwrites the model from
 * it. `setWarmStart(NULL)` refreshes that cache from the model, which is what
 * makes the file's basis survive into the solve. It is also what makes
 * `si.getWarmStart()` -- which is where CglGomory reads the basis at :71 --
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
 * candidate order and for Gomory is a different problem.
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

  // Before the generator: integrality gates candidates, the coefficient formula
  // and the rhs relaxation.
  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  if (!fileExists(bas)) {
    fprintf(stderr, "ERROR: %s: no .bas. Gomory cuts are rows of the tableau at a "
                    "specific basis, so there is nothing to replay without it\n",
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
  // installs the stale cache over it -- and getWarmStart() would hand CglGomory
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
  // same one, which for Gomory is a different algorithm input -- so the CSV
  // carries warmStartIters and a sweep should filter on it.
  if (f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured basis "
                    "did not survive, so these are cuts from a DIFFERENT basis and "
                    "not comparable across runs\n",
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
    "Passing any one of those files also works. The .bas is REQUIRED: a Gomory cut\n"
    "is a row of the simplex tableau at a specific basis.\n"
    "\n"
    "Call shape (defaults reproduce CBC's root call -- see the file comment):\n"
    "  --rounds=N          Gomory rounds, LP re-solved between them (default 1;\n"
    "                      CBC calls Gomory once per pass with a fresh generator)\n"
    "  --pass=N            info.pass (default 0). NOT cosmetic: pass 0 takes\n"
    "                      limitAtRoot (1000) rather than limit (50) AND makes the\n"
    "                      accuracy test at CglGomory.cpp:1252 absolute rather than\n"
    "                      relative -- a different acceptance rule.\n"
    "  --in-tree           info.inTree = true, which switches `away` from\n"
    "                      min(away,awayAtRoot) to away and the limit to limit_\n"
    "  --options=N         info.options bitmask (default 0, which is what\n"
    "                      CbcCutGenerator assigns at the root). Bit 16 makes cuts\n"
    "                      global and sets testFixed=-1; bit 256 is doSorted, which\n"
    "                      turns on the element budget, the candidate sort AND the\n"
    "                      secondaryCuts flush -- so a nonzero value here changes\n"
    "                      the output, not just the speed.\n"
    "\n"
    "CBC's knobs (defaults are CBC's, from CbcSolverCutSetup.cpp:132-237, NOT\n"
    "CglGomory's constructor defaults, which differ):\n"
    "  --limit-at-root=N   cut length cap at the root (default 1000; CBC uses 2000\n"
    "                      above 5000 columns -- pass it explicitly for those).\n"
    "                      NOTE 0 does not disable: :852 turns 0 into numberColumns\n"
    "                      and :876 adds numberRows, i.e. 0 is the LOOSEST setting.\n"
    "  --limit=N           cut length cap in the tree (default 50)\n"
    "  --away-at-root=F    root fractionality threshold (default 0.005, CBC's\n"
    "                      MORE_CUTS value). Usable as a control flag at 0.5,\n"
    "                      where on most instances almost no candidate is admitted,\n"
    "                      so time(0.5) approximates the irreducible per-call floor\n"
    "                      (factorization + rowActivity + rowType) and\n"
    "                      time(0.005) - time(0.5) prices the candidate loop.\n"
    "                      CHECK rowCuts_n BEFORE BELIEVING IT: the test is\n"
    "                      away < f < 1-away, so a basic integer sitting at exactly\n"
    "                      0.5 passes at ANY away<=0.5. On chromaticindex1024-7 it\n"
    "                      still admitted 21314 of 21729 -- graph-colouring LPs are\n"
    "                      half-integral, so the flag prices nothing there. Use\n"
    "                      --no-bound to time separation on those.\n"
    "  --away=F            in-tree threshold (default 0.05)\n"
    "  --gomory-type=N     0 normal (default). 1/2 need passInOriginalSolver to do\n"
    "                      anything, so use --orig-solver with them.\n"
    "  --orig-solver       passInOriginalSolver(clone of the fixture solver), which\n"
    "                      is what makes --gomory-type live. CBC's root does NOT do\n"
    "                      this; it is here to price the path.\n"
    "  --alt-factorization use CoinOslFactorization instead of CoinFactorization\n"
    "  --condition-multiplier=F   default 1e-18 (CBC's), feeds the rhs relaxation\n"
    "  --largest-factor-multiplier=F  default 1e-13 (CBC's), likewise\n"
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
  "                      the re-solve can dwarf what is being measured: on\n"
  "                      chromaticindex1024-7 sepTime is 23.6s and resolveTime is\n"
  "                      1010s, i.e. 96%% of the wall, so a serial min-of-3 on\n"
  "                      sepTime costs ~17 min a rep without it and ~25 s with it.\n"
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

  // scoreColCut must count only genuine tightenings, and must separate "smaller
  // relaxation" from "current point now infeasible". This matters more here than
  // in probing-bench: Gomory reaches OsiColCut through the singleton branch at
  // CglGomory.cpp:1731-1793, so on some fixtures the column cuts are the entire
  // output and a miscount would read as "generated nothing".
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

  // The setters must round-trip, because CglGomory silently ignores out-of-range
  // values: setAway/setAwayAtRoot keep the old value unless 0 < v <= 0.5
  // (CglGomory.cpp:1964-1983) and setLimit/setLimitAtRoot unless v >= 0. A
  // rejected --away would otherwise measure the default while the CSV reported
  // what was asked for, which is exactly the shape of a false result.
  {
    CglGomory g;
    g.setAwayAtRoot(0.005);
    g.setAway(0.05);
    g.setLimitAtRoot(1000);
    g.setLimit(50);
    if (fabs(g.getAwayAtRoot() - 0.005) > 1.0e-15
      || fabs(g.getAway() - 0.05) > 1.0e-15
      || g.getLimitAtRoot() != 1000 || g.getLimit() != 50) {
      fprintf(stderr, "FAIL CBC's root knobs do not round-trip: awayAtRoot=%g "
                      "away=%g limitAtRoot=%d limit=%d\n",
        g.getAwayAtRoot(), g.getAway(), g.getLimitAtRoot(), g.getLimit());
      ++failures;
    }
    // The control flag must be a value the setter accepts, or --away-at-root=0.5
    // would silently leave 0.005 in place and "the floor" would be a full run.
    g.setAwayAtRoot(0.5);
    if (fabs(g.getAwayAtRoot() - 0.5) > 1.0e-15) {
      fprintf(stderr, "FAIL setAwayAtRoot(0.5) rejected -- the control flag is dead, "
                      "got %g\n",
        g.getAwayAtRoot());
      ++failures;
    }
    // And a value outside the accepted range must be REFUSED, not clamped: the
    // CSV prints what was requested, so silent clamping would misreport the run.
    g.setAwayAtRoot(0.9);
    if (fabs(g.getAwayAtRoot() - 0.5) > 1.0e-15) {
      fprintf(stderr, "FAIL setAwayAtRoot(0.9) changed the value to %g; the bench "
                      "assumes out-of-range is ignored\n",
        g.getAwayAtRoot());
      ++failures;
    }
  }

  // needsOptimalBasis() is the premise of this whole bench -- it is why the .bas
  // is required rather than optional, and why a nonzero warmStartIters
  // invalidates a comparison. Pin it: if it ever returns false, that reasoning
  // needs revisiting rather than silently persisting in the comments.
  {
    CglGomory g;
    if (!g.needsOptimalBasis()) {
      fprintf(stderr, "FAIL needsOptimalBasis() is now false; this bench's insistence "
                      "on the captured basis needs rethinking\n");
      ++failures;
    }
  }

  // doSorted, and therefore the secondaryCuts flush, must still be bit 256 of
  // info.options -- the file comment claims those cuts are computed and discarded
  // at CBC's root, which is only true while options==0 fails this test.
  {
    const char *src = "../../Cgl/src/CglGomory/CglGomory.cpp";
    FILE *fp = fopen(src, "r");
    if (!fp) {
      printf("  (skipped doSorted bit check: %s not readable from cwd)\n", src);
    } else {
      char line[512];
      bool found = false;
      while (fgets(line, sizeof(line), fp)) {
        int bit = 0;
        if (sscanf(line, " bool doSorted = (infoOptions&%d)", &bit) == 1) {
          found = true;
          if (bit != 256) {
            fprintf(stderr, "FAIL doSorted is now bit %d, not 256; --options docs and "
                            "the secondaryCuts reasoning are stale\n",
              bit);
            ++failures;
          }
          break;
        }
      }
      fclose(fp);
      if (!found) {
        fprintf(stderr, "FAIL no 'bool doSorted = (infoOptions&N)' in %s\n", src);
        ++failures;
      }
    }
  }

  printf("self-test: %d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}

static const char *CSV_HEADER
  = "name,limitAtRoot,limit,awayAtRoot,away,gomoryType,origSolver,altFactorization,"
    "condMultiplier,largestFactorMultiplier,pass,inTree,options,rounds,"
    "rowCuts_n,colCuts_n,boundsTightened,rowsAdded,totalViol,maxViol,colCutSlack,"
    "colCutViol,avgCutLen,maxCutLen,sepTime,warmStartTime,warmStartIters,"
    "resolveTime,resolveIters,objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "rows,cols,intCols,elements,gomoryCandidates,basicIntegers,nonbasicMovable,"
    "nonbasicMovableNz,denseLoopWork,activeSlackRows,integerSlackRows,restoredInt,"
    "padDropped,rowCutsPerRound,colCutsPerRound,objImprovePerRound";

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
  // One round by default: CBC constructs a fresh CglGomory per cut pass, so a
  // single round is the call being optimized.
  int maxRounds = 1;

  // CBC's values (CbcSolverCutSetup.cpp:132-237), not CglGomory's constructor
  // defaults. Note limitAtRoot: CBC uses 2000 when the *preprocessed* model has
  // more than 5000 columns, which the fixture's .meta records as cbcLimitAtRoot;
  // the default here is the common 1000 and a sweep over large models should pass
  // the meta's value explicitly.
  int limitAtRoot = 1000;
  int limit = 50;
  double awayAtRoot = 0.005;
  double away = 0.05;
  int gomoryType = 0;
  bool origSolver = false;
  bool altFactorization = false;
  double condMultiplier = 1.0e-18;
  double largestFactorMultiplier = 1.0e-13;

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
    } else if (strcmp(a, "--orig-solver") == 0) {
      origSolver = true;
    } else if (strcmp(a, "--alt-factorization") == 0) {
      altFactorization = true;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--pass=", 7) == 0) {
      pass = atoi(a + 7);
    } else if (strncmp(a, "--options=", 10) == 0) {
      options = atoi(a + 10);
    } else if (strncmp(a, "--limit-at-root=", 16) == 0) {
      limitAtRoot = atoi(a + 16);
    } else if (strncmp(a, "--limit=", 8) == 0) {
      limit = atoi(a + 8);
    } else if (strncmp(a, "--away-at-root=", 15) == 0) {
      awayAtRoot = atof(a + 15);
    } else if (strncmp(a, "--away=", 7) == 0) {
      away = atof(a + 7);
    } else if (strncmp(a, "--gomory-type=", 14) == 0) {
      gomoryType = atoi(a + 14);
    } else if (strncmp(a, "--condition-multiplier=", 23) == 0) {
      condMultiplier = atof(a + 23);
    } else if (strncmp(a, "--largest-factor-multiplier=", 28) == 0) {
      largestFactorMultiplier = atof(a + 28);
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

  Fixture f;
  if (!loadFixture(f, stem, quiet))
    return 1;

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();
  const int nElements0 = f.si.getNumElements();

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
    // A fresh generator per round, matching CBC, which constructs one per cut
    // pass. It also removes any chance of state carrying between rounds and
    // confounding a sweep -- the trap that made the clique bench do the same.
    // For Gomory the carried state would be numberTimesStalled_ (dead, both
    // readers commented out) and originalSolver_, but relying on that staying
    // true is exactly the assumption a fresh object removes.
    CglGomory gomory;
    gomory.setLimitAtRoot(limitAtRoot);
    gomory.setLimit(limit);
    gomory.setAwayAtRoot(awayAtRoot);
    gomory.setAway(away);
    gomory.setGomoryType(gomoryType);
    gomory.setConditionNumberMultiplier(condMultiplier);
    gomory.setLargestFactorMultiplier(largestFactorMultiplier);
    if (altFactorization)
      gomory.useAlternativeFactorization(true);
    // Clones internally, so this is a per-round cost and deliberately inside the
    // loop rather than hoisted -- CBC would pay it per pass too, if it did this
    // at all (it does not).
    if (origSolver)
      gomory.passInOriginalSolver(&f.si);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve. The bounds are kept for the
    // same reason -- a column cut is scored against the bounds it tightened.
    const int nc = f.si.getNumCols();
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + nc);
    const std::vector< double > loRound(f.si.getColLower(), f.si.getColLower() + nc);
    const std::vector< double > upRound(f.si.getColUpper(), f.si.getColUpper() + nc);

    OsiCuts cs;

    // The Osi-level entry, which is what CbcCutGenerator calls (it has no
    // Gomory special case, unlike Probing). With gomoryType_==0 and no original
    // solver this is a pass-through to the workhorse overload, so calling the
    // public entry costs nothing and keeps the intVar[] construction and the
    // getWarmStart() call inside the timed region where CBC has them.
    CglTreeInfo info;
    info.level = 0;
    info.pass = pass;
    info.formulation_rows = nRows0;
    info.inTree = inTree;
    info.options = options;

    GOM_PROF_RESET();
    const double t0 = wallClock();
    gomory.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;
    {
      // Tagged with the fixture and round so a multi-fixture serial sweep can be
      // reduced with awk without losing which row each stage belongs to.
      char tag[256];
      snprintf(tag, sizeof(tag), "%s r%d", baseName(stem).c_str(), round);
      GOM_PROF_PRINT(tag);
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
      // Reported alongside the mean because `limit` caps exactly this, so a
      // sweep over --limit-at-root wants to see whether the cap was reached at
      // all before concluding the knob did nothing.
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
    // nothing next to the re-solve. Only the re-solve is skipped, so isProvenOptimal
    // still holds from the warm-started LP and the bound below reads as no movement,
    // which is what the flag documents.
    if (!noBound) {
      const double t1 = wallClock();
      f.si.resolve();
      totalResolveTime += wallClock() - t1;
      totalResolveIters += f.si.getIterationCount();
    }

    // A round whose LP does not reach optimality contributes 0 rather than a
    // bound taken from an unsolved LP. Note infeasibility is a real outcome here
    // and not a failure: the singleton branch emits a 1<=0 row cut precisely to
    // say the node is dead (CglGomory.cpp:1765-1773).
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

  // Unconditional, and not silenced by --quiet: the row about to be printed carries
  // objImprove=0 and resolveTime=0, which is indistinguishable from a change that
  // destroyed the bound. Anyone reading such a row out of a log needs to know why.
  if (noBound)
    fprintf(stderr, "[gomory-bench] --no-bound: objEnd/objImprove/boundMoved and "
                    "resolveTime/resolveIters are NOT measured (reported as no "
                    "movement). Timing of sepTime only.\n");

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%g,%g,%d,%d,%d,%g,%g,%d,%d,%d,%d,"
         "%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g,%.3f,%d,%.6f,%.6f,%d,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%.0f,%ld,%ld,%d,%d,"
         "%s,%s,%s\n",
    baseName(stem).c_str(), limitAtRoot, limit, awayAtRoot, away, gomoryType,
    (int)origSolver, (int)altFactorization, condMultiplier,
    largestFactorMultiplier, pass - round, (int)inTree, options, round,
    totalRowCuts, totalColCuts, totalTightened, f.si.getNumRows() - nRows0,
    totalViol, maxViol, totalColSlack, totalColViol,
    totalRowCuts ? (double)totalCutLen / totalRowCuts : 0.0, maxCutLen,
    totalSepTime, f.warmStartTime, f.warmStartIters,
    totalResolveTime, totalResolveIters, objStart, objEnd, objImprove,
    objImproveRel, (int)boundMoved(objStart, objImprove),
    nRows0, f.si.getNumCols(), intCols, nElements0,
    metaInt(meta, "gomoryCandidates", -1), metaInt(meta, "basicIntegers", -1),
    metaInt(meta, "nonbasicMovable", -1), metaInt(meta, "nonbasicMovableNz", -1),
    metaNum(meta, "denseLoopWork", -1.0),
    metaInt(meta, "activeSlackRows", -1), metaInt(meta, "integerSlackRows", -1),
    f.restoredColTypes, (int)f.paddedRowDropped,
    rowCutsPerRound.c_str(), colCutsPerRound.c_str(), objImprovePerRound.c_str());

  return 0;
}
