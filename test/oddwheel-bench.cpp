/**
 * Stand-alone benchmark for odd-wheel (odd-cycle) separation.
 *
 * @file oddwheel-bench.cpp
 * @brief replay a CglOddWheel call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcClqFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, serialized conflict graph -- warm-starts the LP, and runs N
 * rounds of CglOddWheel. The fixtures were captured for clique separation, and
 * they serve here unchanged because odd-wheel separation consumes exactly the
 * same three things: the conflict graph, the LP solution and the reduced costs.
 * Reaching that state inside CBC costs a full solve, which is what this exists
 * to avoid paying on every iteration of an optimization loop.
 *
 * Figures of merit, in order of authority -- the same order as bkclique-bench,
 * and for the same reasons:
 *
 *   - **bound improvement on reoptimizing** (`objImprove`). The only figure that
 *     says the cuts tightened the relaxation, which is the point of finding them.
 *   - total **violation** (`totalViol`), the proxy available before a re-solve.
 *   - separation **time**, since separation runs many times per solve.
 *
 * The cut count is reported and is the weakest of the four: more cuts at equal
 * bound movement is strictly worse.
 *
 * ## Why this drives CglOddWheel and not CoinOddWheelSeparator
 *
 * The separator could be constructed directly, and that would be a mistake. A
 * standalone probe measures whatever structure the probe itself builds, which is
 * not necessarily the one the real call site builds -- for clique separation a
 * probe of this shape reported 0.95s against ~19ms actually spent, because it
 * built the full conflict graph where the real path uses a small induced
 * subgraph. So the generator is driven through its public entry point, and the
 * separator's internals are read back through CglOddWheel::stats(), which
 * carries CoinOddWheelSeparator::Stats from the last call.
 *
 * ## The one flag that matters most
 *
 * `--ext-method=0` disables the wheel-center lifting entirely
 * (`searchOddWheels()` gates the lifting loop on `extMethod_ > 0`), so
 * `--ext-method=2` minus `--ext-method=0` isolates the cost of lifting exactly.
 * The analogous control is what located the real bottleneck in clique
 * separation in a single sweep -- there it was extension, not the search.
 *
 * ## Reproducibility
 *
 * `--max-seconds` defaults to 0, i.e. no limit, deliberately. CbcCutGenerator
 * hands the generator whatever wall clock remains, and under such a budget the
 * separator aborts partway through graph preparation on a large graph -- so a
 * nonzero limit makes the result depend on machine load and the run stops being
 * reproducible. With the limit off, the counters are a function of the fixture
 * alone.
 *
 * There is no `--min-viol`: the threshold is a compile-time constant in
 * CoinOddWheelSeparator.cpp (ODDWHEEL_SEP_DEF_MIN_VIOL) with no setter, and
 * inventing one just to have a flag would change the code under measurement.
 *
 * The graph is loaded rather than rebuilt by default; a rebuild does not
 * reproduce the graph CBC separated on (measured: trdtaunimep 190900 direct
 * conflicts rebuilt against 191028 used). --rebuild-cgraph selects the old
 * behaviour so the two can be compared directly.
 *
 * Usage:
 *   oddwheel-bench <stem> [options]          (stem = dir/name.tag)
 *   oddwheel-bench --self-test <stem>
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.cgraph, <stem>.bas.
 * A full path to any one of those files also works; the suffix is stripped.
 */

#include "CglOddWheel.hpp"
#include "ClpSimplex.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>

static int g_failures = 0;

static void check(bool ok, const char *msg)
{
  if (!ok) {
    printf("FAIL: %s\n", msg);
    ++g_failures;
  }
}

/// Wall-clock, not CPU: separation cost as a user would feel it.
static double wallClock()
{
  return CoinGetTimeOfDay();
}

static bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

/**
 * Reduce any of the fixture's file names to the shared stem, so a caller can
 * pass whichever one tab-completion produced.
 */
static std::string fixtureStem(const char *arg)
{
  std::string s(arg);
  static const char *suffixes[]
    = { ".mps.gz", ".mps", ".cgraph", ".bas", ".sol", ".bas.status" };
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
 * Read one integer key out of the fixture's `.meta`.
 * Returns `dflt` when the file or the key is absent, so a fixture written before
 * `.meta` carried that key still loads.
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
 * which writeMps would drop those columns and shift every index after them --
 * see cbcClqFixtureWriteMps. The row is redundant, so leaving it in place would
 * not change the LP, but it would change the row count the separator sees and
 * so the cut and row-count columns; drop it and the model is the captured one.
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
 * MPS cannot express "fixed and integer": CoinMpsIO::writeMps conveys
 * integrality only through the bound type (BV/UI/LI/MI+UI) and a column with
 * lb == ub takes the " FX " branch, which has no integer form -- so every integer
 * column CBC had fixed by bound tightening reads back continuous.
 *
 * CglOddWheel itself reads only the conflict graph, so integrality does not gate
 * it the way it gates the clique-row detector. It is restored anyway, because
 * refreshSolver() recomputes column types and because a fixture that loads
 * differently here than in bkclique-bench would not be the same experiment.
 *
 * Returns the number of columns re-marked, or -1 when no usable sidecar was
 * found. A sidecar whose column count disagrees with the loaded model is refused
 * rather than partly applied.
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
  // it has to be recomputed after this or consumers would see the pre-restore view.
  si.getColType(true);
  return restored;
}

/// Everything the fixture determines, plus how long loading it took.
struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  /// Pivots the warm start needed. Zero is the expected value and the check that
  /// the captured basis actually took; see loadFixture below.
  int warmStartIters = 0;
  double cgraphTime = 0.0;
  bool paddedRowDropped = false;
  int restoredColTypes = -1;
  bool ok = false;
};

/**
 * Load the fixture and warm-start to the captured optimum.
 *
 * Two things are needed to actually *land* on the captured vertex, and getting
 * either wrong looks like success while silently changing the experiment -- the
 * odd cycles found are those violated by the current fractional point, so a
 * different optimal vertex means a different set of cuts.
 *
 * First, `readBasis` writes into ClpSimplex's own `status_` array, but OsiClp
 * caches a separate `CoinWarmStartBasis basis_` and both entry points overwrite
 * the model from it -- `resolve()` at OsiClpSolverInterface.cpp:1199, and
 * `initialSolve()` by presolving from scratch. `setWarmStart(NULL)` refreshes
 * `basis_` from the model (`basis_ = getBasis(modelPtr_)`, applying the slack
 * flip), which is what makes the file's basis survive into the solve.
 *
 * Second, `resolve()` rather than `initialSolve()`: presolve discards the basis.
 * With both in place an already-optimal fixture costs 0 iterations, which is the
 * cheap self-check that the warm start worked. `setPerturbation(50)` matches the
 * generator; perturbation left on moves the vertex even from a correct basis
 * (decomp2: 846 columns moved versus 803).
 */
static bool loadFixture(Fixture &f, const std::string &stem, bool rebuildCgraph,
  bool quiet)
{
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";
  const std::string cgr = stem + ".cgraph";

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

  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  bool haveBasis = false;
  if (fileExists(bas)) {
    if (lp->readBasis(bas.c_str()) < 0) {
      fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
    } else {
      haveBasis = true;
      // Push the model's freshly-read status into OsiClp's cached basis_, or the
      // solve below installs the stale cache over it.
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
  // A correct warm start from an optimal basis costs no pivots. Anything else means
  // the fixture landed on a different vertex, so say so rather than quietly
  // measuring a different LP.
  if (haveBasis && f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured basis "
                    "did not survive, so this is a different vertex\n",
      baseName(stem).c_str(), f.warmStartIters);
  }

  const double t1 = wallClock();
  if (rebuildCgraph) {
    f.si.checkCGraph(NULL);
  } else {
    if (!fileExists(cgr)) {
      fprintf(stderr, "ERROR: no conflict graph %s (use --rebuild-cgraph to build one)\n",
        cgr.c_str());
      return false;
    }
    CoinStaticConflictGraph *cg = CoinStaticConflictGraph::load(cgr.c_str());
    if (!cg) {
      fprintf(stderr, "ERROR: failed to load conflict graph %s\n", cgr.c_str());
      return false;
    }
    if (cg->size() != (size_t)f.si.getNumCols() * 2) {
      // CglOddWheel calls exit() on this, so catch it here with a message that
      // says which fixture is inconsistent.
      fprintf(stderr, "ERROR: graph/model mismatch for %s: graph %lu nodes, model %d cols\n",
        stem.c_str(), (unsigned long)cg->size(), f.si.getNumCols());
      delete cg;
      return false;
    }
    f.si.setCGraph(cg);
  }
  f.cgraphTime = wallClock() - t1;

  f.ok = true;
  return true;
}

/**
 * Round-trip the graph through save()/load() and prove the reconstruction is
 * identical, field by field. "Close" is not good enough here: the whole point of
 * serializing is that a fixture reproduces what the separator saw exactly,
 * including the deliberately *approximate* degrees, which a rebuild would
 * recompute differently -- and degree gates which nodes are active at all
 * (fillActiveColumns skips degree < 2) as well as which hole node is picked as
 * the wheel-center seed.
 */
static int selfTest(const std::string &stem)
{
  const std::string cgr = stem + ".cgraph";
  if (!fileExists(cgr)) {
    fprintf(stderr, "ERROR: no conflict graph %s\n", cgr.c_str());
    return 1;
  }

  CoinStaticConflictGraph *a = CoinStaticConflictGraph::load(cgr.c_str());
  if (!a) {
    printf("FAIL: could not load %s\n", cgr.c_str());
    return 1;
  }

  const std::string t1 = "/tmp/oddwheel-selftest-1.cgraph";
  const std::string t2 = "/tmp/oddwheel-selftest-2.cgraph";
  check(a->save(t1.c_str()) == 0, "save of loaded graph");

  CoinStaticConflictGraph *b = CoinStaticConflictGraph::load(t1.c_str());
  if (!b) {
    printf("FAIL: could not re-load saved graph\n");
    delete a;
    return 1;
  }
  check(b->save(t2.c_str()) == 0, "re-save of re-loaded graph");

  // Byte identity of save(load(save(x))) vs save(x): catches any field that is
  // written but silently not read back, which a getter-by-getter comparison
  // would miss if the getter does not expose it.
  {
    FILE *f1 = fopen(t1.c_str(), "rb");
    FILE *f2 = fopen(t2.c_str(), "rb");
    bool same = f1 && f2;
    long n = 0;
    while (same) {
      char b1[65536], b2[65536];
      const size_t r1 = fread(b1, 1, sizeof(b1), f1);
      const size_t r2 = fread(b2, 1, sizeof(b2), f2);
      if (r1 != r2 || memcmp(b1, b2, r1) != 0) {
        same = false;
        break;
      }
      n += (long)r1;
      if (r1 == 0)
        break;
    }
    if (f1)
      fclose(f1);
    if (f2)
      fclose(f2);
    check(same, "save/load/save is byte-identical");
    if (same)
      printf("  round-trip stable over %ld bytes\n", n);
  }

  check(a->size() == b->size(), "size preserved");
  check(a->density() == b->density(), "density preserved");
  check(a->minDegree() == b->minDegree(), "minDegree preserved");
  check(a->maxDegree() == b->maxDegree(), "maxDegree preserved");
  check(a->nCliques() == b->nCliques(), "nCliques preserved");
  check(a->nTotalDirectConflicts() == b->nTotalDirectConflicts(),
    "total direct conflicts preserved");
  check(a->nTotalCliqueElements() == b->nTotalCliqueElements(),
    "total clique elements preserved");

  if (a->size() == b->size()) {
    size_t badDeg = 0, badMDeg = 0, badConf = 0, badNodeClq = 0;
    for (size_t i = 0; i < a->size(); ++i) {
      if (a->degree(i) != b->degree(i))
        ++badDeg;
      if (a->modifiedDegree(i) != b->modifiedDegree(i))
        ++badMDeg;
      if (a->nNodeCliques(i) != b->nNodeCliques(i)) {
        ++badNodeClq;
        continue;
      }
      for (size_t k = 0; k < a->nNodeCliques(i); ++k)
        if (a->nodeCliques(i)[k] != b->nodeCliques(i)[k])
          ++badNodeClq;
      if (a->nDirectConflicts(i) != b->nDirectConflicts(i)) {
        ++badConf;
        continue;
      }
      for (size_t k = 0; k < a->nDirectConflicts(i); ++k)
        if (a->directConflicts(i)[k] != b->directConflicts(i)[k])
          ++badConf;
    }
    check(badDeg == 0, "per-node degree preserved");
    check(badMDeg == 0, "per-node modified degree preserved");
    check(badConf == 0, "per-node direct conflicts preserved");
    check(badNodeClq == 0, "per-node clique index rebuilt correctly");
  }

  if (a->nCliques() == b->nCliques()) {
    size_t badClq = 0;
    for (size_t c = 0; c < a->nCliques(); ++c) {
      if (a->cliqueSize(c) != b->cliqueSize(c)) {
        ++badClq;
        continue;
      }
      for (size_t k = 0; k < a->cliqueSize(c); ++k)
        if (a->cliqueElements(c)[k] != b->cliqueElements(c)[k])
          ++badClq;
    }
    check(badClq == 0, "clique elements preserved");
  }

  printf("%s: size=%lu directConf=%lu cliques=%lu clqElems=%lu density=%.10f\n",
    baseName(stem).c_str(), (unsigned long)a->size(),
    (unsigned long)a->nTotalDirectConflicts(), (unsigned long)a->nCliques(),
    (unsigned long)a->nTotalCliqueElements(), a->density());

  delete a;
  delete b;

  printf(g_failures ? "self-test: %d FAILURE(S)\n" : "self-test: all checks passed\n",
    g_failures);
  return g_failures ? 1 : 0;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "       %s --self-test <fixture-stem>\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.cgraph/.bas are appended.\n"
    "Passing any one of those files also works.\n"
    "\n"
    "Options:\n"
    "  --rounds=N          separation rounds, LP re-solved between them (default 4)\n"
    "  --ext-method=N      wheel-center lifting: 0 = none, 1 = one variable,\n"
    "                      2 = a clique (default 2, CBC's). 0 isolates the cost\n"
    "                      of lifting, since searchOddWheels() skips the whole\n"
    "                      lifting loop when it is 0.\n"
    "  --max-seconds=F     separator wall-clock budget (default 0 = none; a\n"
    "                      nonzero value makes the run load-dependent)\n"
    "  --rebuild-cgraph    rebuild the graph from the matrix instead of loading\n"
    "                      the captured one (not faithful; for comparison only)\n"
    "  --stage-times       also print a human-readable stage breakdown to stderr\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --quiet             suppress warnings\n",
    prog, prog);
}

/**
 * Did the bound actually move, or is this floating-point noise?
 *
 * objImprove is a difference of two LP objectives, so it inherits their absolute
 * error, which scales with their magnitude -- an LP around 1e6 carries roughly
 * 1e-10 of slop. Reading such a difference against zero therefore reports an
 * improvement on every fixture: measured, brazil3's clqstr-before pass gives
 * 2.7e-12 over 29317 iterations, which is nothing at all, while mzzv11 gives
 * 171.49, which is real. A relative test separates the two, with an absolute
 * floor so that an objective near zero does not make the relative test
 * hypersensitive.
 */
static bool boundMoved(double objStart, double objImprove)
{
  const double scale = fabs(objStart) > 1.0 ? fabs(objStart) : 1.0;
  return objImprove > 1.0e-9 * scale && objImprove > 1.0e-9;
}

static const char *CSV_HEADER
  = "name,extMethod,rounds,rowsAdded,totalCuts,totalViol,maxViol,avgCutLen,"
    "sepTime,warmStartTime,warmStartIters,cgraphTime,resolveTime,resolveIters,"
    "objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "cgNodes,cgDirectConf,cgCliques,cgDensity,restoredInt,"
    "icaCount,arcs,spFindCalls,oddHoles,ohShort,ohRepeated,ohNotViol,ohDuplicate,"
    "wheelCenters,wcElements,avgWcSize,cutsBeforePool,cutsDupIdx,cutsZeroCoefs,"
    "cutsEmpty,cutsAfterPool,"
    "timeLimitHit,tSetup,tSeparator,tActive,tPrepArcs,tPrepRev,tPrepSpf,tSearch,"
    "tWheelCenter,tCutPool,cutsPerRound,violPerRound,objImprovePerRound";

/**
 * Sums of CglOddWheel::stats() over the rounds.
 *
 * Kept separate from the per-call struct rather than accumulated into it: the
 * generator is rebuilt every round (see the loop) so its own stats are always
 * one call's worth, which is what makes a single round's numbers readable in
 * isolation when --rounds=1.
 */
struct Totals {
  size_t icaCount = 0, arcs = 0, spFindCalls = 0;
  size_t oddHoles = 0, ohShort = 0, ohRepeated = 0, ohNotViol = 0, ohDuplicate = 0;
  size_t wheelCenters = 0, wcElements = 0;
  size_t cutsBeforePool = 0, cutsDupIdx = 0, cutsAfterPool = 0;
  size_t cutsZeroCoefs = 0, cutsEmpty = 0;
  bool timeLimitHit = false;
  double tSetup = 0.0, tSeparator = 0.0, tActive = 0.0, tPrepArcs = 0.0;
  double tPrepRev = 0.0, tPrepSpf = 0.0, tSearch = 0.0, tWheelCenter = 0.0;
  double tCutPool = 0.0;

  void add(const CglOddWheel::Stats &s)
  {
    // icaCount is a size, not a count of events, so the max over the rounds is
    // the meaningful summary: rounds after the first separate a *larger* LP but
    // usually a less fractional one, and summing would report a node count no
    // single call ever saw.
    if (s.sep.activeColumns > icaCount)
      icaCount = s.sep.activeColumns;
    if (s.sep.arcs > arcs)
      arcs = s.sep.arcs;

    spFindCalls += s.sep.spFindCalls;
    oddHoles += s.sep.oddHolesFound;
    ohShort += s.sep.oddHolesShort;
    ohRepeated += s.sep.oddHolesRepeatedNode;
    ohNotViol += s.sep.oddHolesNotViolated;
    ohDuplicate += s.sep.oddHolesDuplicate;
    wheelCenters += s.sep.wheelCenters;
    wcElements += s.sep.wheelCenterElements;
    timeLimitHit = timeLimitHit || s.sep.timeLimitReached;

    cutsBeforePool += s.cutsBeforePool;
    cutsDupIdx += s.cutsDuplicatedIdx;
    cutsZeroCoefs += s.cutsZeroCoefs;
    cutsEmpty += s.cutsEmpty;
    cutsAfterPool += s.cutsAfterPool;

    tSetup += s.tSetup;
    tSeparator += s.tSeparator;
    tActive += s.sep.tActiveColumns;
    tPrepArcs += s.sep.tPrepareArcs;
    tPrepRev += s.sep.tPrepareReverse;
    tPrepSpf += s.sep.tPrepareShortestPath;
    tSearch += s.sep.tSearch;
    tWheelCenter += s.sep.tWheelCenter;
    tCutPool += s.tCutPool;
  }
};

/**
 * The stage breakdown, on stderr so it never contaminates the CSV line.
 *
 * tSearch is reported as one figure and deliberately not split per
 * CoinShortestPath::find() call: two clock reads times icaCount_ calls is
 * hundreds of thousands of reads on a large fixture, which would perturb the
 * very loop being measured. The Dijkstra share comes out by difference instead.
 */
static void printStageTimes(const Totals &t, double totalSepTime)
{
  const double stages = t.tSetup + t.tActive + t.tPrepArcs + t.tPrepRev + t.tPrepSpf
    + t.tSearch + t.tWheelCenter + t.tCutPool;
  const double pct = totalSepTime > 0.0 ? 100.0 / totalSepTime : 0.0;

  fprintf(stderr, "\n  stage                       seconds     %% of sepTime\n");
  fprintf(stderr, "  doubled x/rc setup       %10.4f   %8.2f\n", t.tSetup, t.tSetup * pct);
  fprintf(stderr, "  fillActiveColumns        %10.4f   %8.2f\n", t.tActive, t.tActive * pct);
  fprintf(stderr, "  prepareGraph: conflicts  %10.4f   %8.2f\n", t.tPrepArcs, t.tPrepArcs * pct);
  fprintf(stderr, "  prepareGraph: reverse    %10.4f   %8.2f\n", t.tPrepRev, t.tPrepRev * pct);
  fprintf(stderr, "  prepareGraph: sp ctor    %10.4f   %8.2f\n", t.tPrepSpf, t.tPrepSpf * pct);
  fprintf(stderr, "  odd-hole search          %10.4f   %8.2f\n", t.tSearch, t.tSearch * pct);
  fprintf(stderr, "  wheel-center lifting     %10.4f   %8.2f\n", t.tWheelCenter,
    t.tWheelCenter * pct);
  fprintf(stderr, "  cut pool + insertion     %10.4f   %8.2f\n", t.tCutPool, t.tCutPool * pct);
  fprintf(stderr, "  ------------------------------------------------\n");
  fprintf(stderr, "  accounted                %10.4f   %8.2f\n", stages, stages * pct);
  fprintf(stderr, "  sepTime (wall)           %10.4f   %8.2f\n", totalSepTime, 100.0);
  fprintf(stderr, "\n  shortest-path calls %lu, odd holes kept %lu of %lu paths"
                  " (%lu short, %lu repeated node, %lu not violated, %lu duplicate)\n",
    (unsigned long)t.spFindCalls, (unsigned long)t.oddHoles,
    (unsigned long)t.spFindCalls, (unsigned long)t.ohShort,
    (unsigned long)t.ohRepeated, (unsigned long)t.ohNotViol,
    (unsigned long)t.ohDuplicate);
  fprintf(stderr, "  wheel centers %lu with %lu elements; cuts %lu pooled,"
                  " %lu merged (%lu coefs cancelled, %lu emptied), %lu survived\n\n",
    (unsigned long)t.wheelCenters, (unsigned long)t.wcElements,
    (unsigned long)t.cutsBeforePool, (unsigned long)t.cutsDupIdx,
    (unsigned long)t.cutsZeroCoefs, (unsigned long)t.cutsEmpty,
    (unsigned long)t.cutsAfterPool);
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
  if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
    usage(argv[0]);
    return 0;
  }

  bool doSelfTest = false;
  bool rebuildCgraph = false;
  bool csvHeader = false;
  bool quiet = false;
  bool stageTimes = false;
  int maxRounds = 4;
  size_t extMethod = 2;
  double maxSeconds = 0.0;
  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (strcmp(a, "--self-test") == 0) {
      doSelfTest = true;
    } else if (strcmp(a, "--rebuild-cgraph") == 0) {
      rebuildCgraph = true;
    } else if (strcmp(a, "--csv-header") == 0) {
      csvHeader = true;
    } else if (strcmp(a, "--quiet") == 0) {
      quiet = true;
    } else if (strcmp(a, "--stage-times") == 0) {
      stageTimes = true;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--ext-method=", 13) == 0) {
      extMethod = (size_t)atol(a + 13);
      if (extMethod > 2) {
        fprintf(stderr, "ERROR: --ext-method must be 0, 1 or 2\n");
        return 1;
      }
    } else if (strncmp(a, "--max-seconds=", 14) == 0) {
      maxSeconds = atof(a + 14);
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

  if (doSelfTest)
    return selfTest(stem);

  Fixture f;
  if (!loadFixture(f, stem, rebuildCgraph, quiet))
    return 1;

  const CoinConflictGraph *cg = f.si.getCGraph();
  if (!cg) {
    fprintf(stderr, "ERROR: no conflict graph available for %s\n", stem.c_str());
    return 1;
  }

  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  info.formulation_rows = f.si.getNumRows();
  info.inTree = false;
  info.options = 0;

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();

  int totalCuts = 0;
  size_t totalCutLen = 0;
  double totalViol = 0.0, maxViol = 0.0;
  double totalSepTime = 0.0, totalResolveTime = 0.0;
  int totalResolveIters = 0;
  Totals tot;
  std::string cutsPerRound, violPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution can be
  // reported: a separator whose whole gain arrives in round 1 behaves differently
  // under Cbc's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A fresh generator per round. CglOddWheel keeps no ratchet of its own, but
    // it does keep the scratch arrays and the stats of the last call, and a fresh
    // object makes each round independent of the previous one by construction
    // rather than by inspection.
    CglOddWheel oddWheel(extMethod);
    if (maxSeconds > 0.0)
      oddWheel.setMaxSeconds(maxSeconds);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve.
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + f.si.getNumCols());

    OsiCuts cs;
    const double t0 = wallClock();
    oddWheel.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;
    tot.add(oddWheel.stats());

    const int nCuts = cs.sizeRowCuts();
    if (nCuts == 0)
      break;

    double roundViol = 0.0;
    for (int c = 0; c < nCuts; ++c) {
      const OsiRowCut &rc = cs.rowCut(c);
      const double v = rc.violated(xRound.data());
      roundViol += v;
      if (v > maxViol)
        maxViol = v;
      totalCutLen += (size_t)rc.row().getNumElements();
    }
    totalViol += roundViol;
    totalCuts += nCuts;

    char buf[64];
    snprintf(buf, sizeof(buf), "%s%d", round ? "+" : "", nCuts);
    cutsPerRound += buf;
    snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "", roundViol);
    violPerRound += buf;

    f.si.applyCuts(cs);

    const double t1 = wallClock();
    f.si.resolve();
    totalResolveTime += wallClock() - t1;
    totalResolveIters += f.si.getIterationCount();

    // Per-round bound gain. A round whose LP does not reach optimality
    // contributes 0 rather than a bound taken from an unsolved LP.
    const double objRoundEnd
      = f.si.isProvenOptimal() ? f.si.getObjValue() : objRoundStart;
    snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "",
      f.si.getObjSenseInCbc() * (objRoundEnd - objRoundStart));
    objImprovePerRound += buf;
    objRoundStart = objRoundEnd;

    info.pass = round + 1;
  }

  if (cutsPerRound.empty()) {
    cutsPerRound = "0";
    violPerRound = "0";
    objImprovePerRound = "0";
  }

  const double objEnd = f.si.isProvenOptimal() ? f.si.getObjValue() : objStart;
  const double objImprove = f.si.getObjSenseInCbc() * (objEnd - objStart);
  // Relative to the starting bound's magnitude, which is what makes gains
  // comparable across instances whose objectives differ by orders of magnitude.
  const double objImproveRel
    = objImprove / (fabs(objStart) > 1.0 ? fabs(objStart) : 1.0);

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%lu,%d,%d,%d,%.10g,%.10g,%.3f,%.6f,%.6f,%d,%.6f,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,"
         "%lu,%lu,%lu,%.10f,%d,"
         "%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,"
         "%lu,%lu,%.3f,%lu,%lu,%lu,%lu,%lu,"
         "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
         "%.6f,%.6f,%s,%s,%s\n",
    baseName(stem).c_str(), (unsigned long)extMethod, round,
    f.si.getNumRows() - nRows0, totalCuts, totalViol, maxViol,
    totalCuts ? (double)totalCutLen / totalCuts : 0.0,
    totalSepTime, f.warmStartTime, f.warmStartIters, f.cgraphTime,
    totalResolveTime, totalResolveIters, objStart, objEnd, objImprove,
    objImproveRel, (int)boundMoved(objStart, objImprove),
    (unsigned long)cg->size(), (unsigned long)cg->nTotalDirectConflicts(),
    (unsigned long)cg->nCliques(), cg->density(), f.restoredColTypes,
    (unsigned long)tot.icaCount, (unsigned long)tot.arcs,
    (unsigned long)tot.spFindCalls, (unsigned long)tot.oddHoles,
    (unsigned long)tot.ohShort, (unsigned long)tot.ohRepeated,
    (unsigned long)tot.ohNotViol, (unsigned long)tot.ohDuplicate,
    (unsigned long)tot.wheelCenters, (unsigned long)tot.wcElements,
    tot.wheelCenters ? (double)tot.wcElements / tot.wheelCenters : 0.0,
    (unsigned long)tot.cutsBeforePool, (unsigned long)tot.cutsDupIdx,
    (unsigned long)tot.cutsZeroCoefs, (unsigned long)tot.cutsEmpty,
    (unsigned long)tot.cutsAfterPool, (int)tot.timeLimitHit,
    tot.tSetup, tot.tSeparator, tot.tActive, tot.tPrepArcs, tot.tPrepRev,
    tot.tPrepSpf, tot.tSearch, tot.tWheelCenter, tot.tCutPool,
    cutsPerRound.c_str(), violPerRound.c_str(), objImprovePerRound.c_str());

  if (stageTimes)
    printStageTimes(tot, totalSepTime);

  return 0;
}
