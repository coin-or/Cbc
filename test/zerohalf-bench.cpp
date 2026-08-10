/**
 * Stand-alone benchmark for 0-1/2 (ZeroHalf) cut separation.
 *
 * @file zerohalf-bench.cpp
 * @brief replay a CglZeroHalf separation call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcZeroHalfFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, LP solution, column types -- warm-starts the LP, and runs N
 * rounds of CglZeroHalf. Reaching this state inside CBC costs a preprocess plus a
 * root LP, so this exists to make the iterate-measure loop on the separation
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
 *   - the cut count is reported and is the **weakest** of the four: more cuts at
 *     equal bound movement is strictly worse, being more rows for the same
 *     tightening. Never rank by it and never compare it against zero.
 *
 * NO AUXILIARY STRUCTURE IS LOADED, AND THAT IS FAITHFUL. The clique fixtures had
 * to serialize CBC's conflict graph because CBC builds it once and caches a stale
 * copy. CglZeroHalf has no such cache: refreshSolver() rebuilds mr_/mc_/mnz_/
 * mtbeg_/mtcnt_/mtind_/mtval_/vlb_/vub_/mrhs_/msense_ from the solver's row-major
 * matrix, column bounds and column types on *every* call, and generateCuts()
 * re-derives the integer bounds on top of that. So reconstructing from the model
 * is not an approximation of what CBC did -- it is exactly what CBC does.
 *
 * THE CONTROL FLAG. `--max-pairs=N` caps the per-row pair count, which is the
 * knob that disables the expensive stage: basic_separation() weakens every pair
 * of odd entries in a row, an O(cnt^2) loop whose body calls best_weakening() on
 * the remaining cnt-2 variables, so the cost of a row is O(cnt^3). `--max-pairs=0`
 * therefore prices that stage exactly, in the way `--ext-method=0` priced clique
 * extension. `--max-frac=N` is the second existing gate, on the count of
 * fractional entries.
 *
 * Usage:
 *   zerohalf-bench <stem> [options]          (stem = dir/name.tag)
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.bas, <stem>.sol,
 * <stem>.ctype, <stem>.meta. A full path to any one of those also works; the
 * suffix is stripped.
 */

#include "CglZeroHalf.hpp"
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
 * cbcZhFixtureWriteMps. The row is redundant so leaving it would not change the
 * LP, but it would change the row count the separator sees, and hence what
 * refreshSolver() builds.
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
 * For ZeroHalf that changes the experiment rather than merely annotating it.
 * refreshSolver() rejects any row touching a continuous column
 * (`if (vlb_[jColumn]==COIN_INT_MAX) { good=false; break; }`), so each lost marker
 * deletes every row that column appears in, shrinking mr_/mnz_ and changing what
 * gets separated. Skipping this step silently measures a smaller problem.
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
  // it must be recomputed or refreshSolver() would see the pre-restore view.
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
 * either wrong looks like success while silently changing the experiment --
 * ZeroHalf separates against the current fractional point, so a different optimal
 * vertex means a different parity ILP and a different cut set.
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

  // Before separation: integrality is the gate refreshSolver() applies per row.
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
  // measuring a different LP.
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
    "Options:\n"
    "  --rounds=N          separation rounds, LP re-solved between them (default 1;\n"
    "                      CBC calls ZeroHalf once per pass with a fresh generator)\n"
    "  --flags=N           CglZeroHalf flags (default 1 = global cuts, CBC's root\n"
    "                      setting; bit 0 clear makes it re-derive integer bounds)\n"
    "  --sparse-threshold=N  active-node count above which the separation graph is\n"
    "                      built sparse (default 8000, CBC's; 0 = always sparse)\n"
    "  --max-pairs=N       skip any row whose odd-entry pair count exceeds N\n"
    "                      (default -1 = no skipping). THE CONTROL FLAG: 0 disables\n"
    "                      the O(cnt^3) pair-weakening stage entirely, so\n"
    "                      time(-1) - time(0) prices that stage exactly.\n"
    "  --max-frac=N        skip any row with more than N fractional entries\n"
    "                      (default -1 = no skipping)\n"
    "  --max-seconds=F     generator time budget (default 0 = none; leave at 0,\n"
    "                      a wall-clock budget makes results irreproducible)\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
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

static const char *CSV_HEADER
  = "name,flags,sparseThreshold,maxPairs,maxFrac,rounds,rowsAdded,totalCuts,"
    "totalViol,maxViol,avgCutLen,sepTime,refreshTime,warmStartTime,warmStartIters,"
    "resolveTime,resolveIters,objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "rows,cols,intCols,zhRows,zhNonzeros,restoredInt,"
    "cutsPerRound,violPerRound,objImprovePerRound";

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

  bool csvHeader = false;
  bool quiet = false;
  // One round by default, unlike the clique bench: CBC constructs a fresh
  // CglZeroHalf per cut pass, so a single round is the call being optimized.
  int maxRounds = 1;
  int flags = 1;
  int sparseThreshold = -1; // -1 = leave the generator's own default
  int maxPairs = -1;
  int maxFrac = -1;
  double maxSeconds = 0.0;
  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (strcmp(a, "--csv-header") == 0) {
      csvHeader = true;
    } else if (strcmp(a, "--quiet") == 0) {
      quiet = true;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--flags=", 8) == 0) {
      flags = atoi(a + 8);
    } else if (strncmp(a, "--sparse-threshold=", 19) == 0) {
      sparseThreshold = atoi(a + 19);
    } else if (strncmp(a, "--max-pairs=", 12) == 0) {
      maxPairs = atoi(a + 12);
    } else if (strncmp(a, "--max-frac=", 11) == 0) {
      maxFrac = atoi(a + 11);
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

  Fixture f;
  if (!loadFixture(f, stem, quiet))
    return 1;

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
  double totalSepTime = 0.0, totalResolveTime = 0.0, totalRefreshTime = 0.0;
  int totalResolveIters = 0;
  std::string cutsPerRound, violPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution is
  // reported separately: a generator whose whole gain arrives in round 1 behaves
  // differently under CBC's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A fresh generator per round, matching CBC, which constructs one per cut
    // pass. It also removes any chance of state carrying between rounds and
    // confounding a sweep -- the trap that made the clique bench do the same.
    CglZeroHalf zeroHalf;
    zeroHalf.setFlags(flags);
    if (sparseThreshold >= 0)
      zeroHalf.setSepGraphSparseThreshold(sparseThreshold);
    zeroHalf.setRowMaxPairCount(maxPairs);
    zeroHalf.setRowMaxFractionalCount(maxFrac);
    // maxSeconds left at 0 by default: Cgl012Cut shrinks its own work when given
    // a wall-clock budget, which would make timings self-referential.
    if (maxSeconds > 0.0)
      zeroHalf.setMaxSeconds(maxSeconds);

    // MANDATORY, and silent if forgotten. generateCuts() is wrapped in
    // `if (mnz_)`, and mnz_ is only ever set by refreshSolver() -- which CBC
    // calls from CbcCutGenerator::generateCuts, not CglZeroHalf itself. Omit
    // this and every fixture reports 0 cuts in 0.000s, which reads as "ZeroHalf
    // has nothing to do here" rather than "the benchmark never ran it".
    // Timed separately from separation: this is the parity-ILP construction,
    // which CBC also pays per call, but it is not the algorithm under study.
    const double tr = wallClock();
    zeroHalf.refreshSolver(&f.si);
    totalRefreshTime += wallClock() - tr;

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve.
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + f.si.getNumCols());

    OsiCuts cs;
    const double t0 = wallClock();
    zeroHalf.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;

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

    // A round whose LP does not reach optimality contributes 0 rather than a
    // bound taken from an unsolved LP.
    const double objRoundEnd
      = f.si.isProvenOptimal() ? f.si.getObjValue() : objRoundStart;
    snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "",
      f.si.getObjSense() * (objRoundEnd - objRoundStart));
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

  printf("%s,%d,%d,%d,%d,%d,%d,%d,%.10g,%.10g,%.3f,%.6f,%.6f,%.6f,%d,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,%d,%d,%d,%ld,%ld,%d,%s,%s,%s\n",
    baseName(stem).c_str(), flags, sparseThreshold, maxPairs, maxFrac,
    round, f.si.getNumRows() - nRows0, totalCuts, totalViol, maxViol,
    totalCuts ? (double)totalCutLen / totalCuts : 0.0,
    totalSepTime, totalRefreshTime, f.warmStartTime, f.warmStartIters,
    totalResolveTime, totalResolveIters, objStart, objEnd, objImprove,
    objImproveRel, (int)boundMoved(objStart, objImprove),
    nRows0, f.si.getNumCols(), intCols,
    metaInt(meta, "zhRows", -1), metaInt(meta, "zhNonzeros", -1),
    f.restoredColTypes,
    cutsPerRound.c_str(), violPerRound.c_str(), objImprovePerRound.c_str());

  return 0;
}
