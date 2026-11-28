/**
 * rins-bench — replay RINS (and, best-effort, VND) from a fixture written by
 * CbcHeuristicRinsFixtureDump.hpp (CBC_DUMP_RINS_FIXTURE), for fast parameter
 * sweeps without paying for a full branch-and-bound run per experiment.
 *
 * @file rins-bench.cpp
 *
 * Loads the current node's LP relaxation (matrix + optimal basis) and the
 * incumbent solution captured at the exact call site inside
 * CbcHeuristicRINS::solution(), rebuilds a CbcModel around it, installs the
 * incumbent as the model's best solution, attaches a REAL CbcHeuristicRINS
 * (the production class, not a reimplementation) configured with the given
 * knobs, and calls heuristic->solution() directly -- exactly the call
 * CbcModel's own search loop makes, just without the surrounding tree. This
 * makes every sweep point cost one sub-MIP solve, not a full CBC run.
 *
 * Also reports, independent of what the heuristic actually did, the
 * "exact-agreement" fix count RINS's own threshold gate is based on (`5*nFix
 * > numberIntegers`) -- this is what lets --sweep rank shallow/close-fix/node
 * combinations by how large a sub-MIP they leave, not just whether a
 * solution was found.
 *
 * See RINS-FIXTURES.md for the fixture format and generation workflow
 * (`gen-rins-fixtures`).
 *
 * Usage:
 *   rins-bench <instance-name> [options]        (looks up the fixture in
 *                                                 --fixture-dir)
 *   rins-bench <stem-or-any-fixture-file> [options]
 *
 * Options:
 *   --fixture-dir=PATH   Fixture directory (default: $CBC_RINS_FIXTURE_DIR or
 *                        ~/instances/mip-sanity-data/rinsFixtures)
 *   --method=rins|vnd    Heuristic to run (default rins)
 *   --shallow=N          shallowDepth_ fixing mode 0-3 (default 0)
 *   --fix-close=X        fixCloseMaxDist_, RINS only (default 0.4)
 *   --nodes=N            Sub-MIP node budget (default 200)
 *   --frac-small=X       fractionSmall_ sub-MIP size gate (default 0.5)
 *   --sec=N              Wall time backstop for the whole replay (default 30)
 *   --probe-report       Additionally run CglProbing on the fixed sub-MIP and
 *                        report how many MORE variables it could fix via
 *                        bound propagation, before any B&B happens
 *   --sweep              Ignore the single-point options above and run a
 *                        built-in grid of shallow x fix-close x nodes,
 *                        printing one row per combination
 *   --log=N              CBC log level for the sub-MIP solve (default 0)
 *   --quiet              Suppress warnings (e.g. missing sidecars)
 *
 * Exit codes: 0 loaded and ran; 1 fixture load failed; 2 usage error.
 */

#include "CbcHeuristicRINS.hpp"
#include "CbcModel.hpp"
#include "ClpSimplex.hpp"
#include "CglProbing.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>

static double wallClock() { return CoinGetTimeOfDay(); }

static bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

static std::string baseName(const std::string &path)
{
  const size_t slash = path.rfind('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

static std::string fixtureStem(const std::string &arg, const std::string &fixtureDir)
{
  static const char *suffixes[] = { ".mps.gz", ".mps", ".bas", ".sol", ".ctype", ".meta",
    ".incumbent", ".dumplog" };
  for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
    const std::string suf(suffixes[i]);
    if (arg.size() > suf.size() && arg.compare(arg.size() - suf.size(), suf.size(), suf) == 0)
      return arg.substr(0, arg.size() - suf.size());
  }
  if (arg.find('/') == std::string::npos)
    return fixtureDir + "/" + arg + ".rins";
  return arg;
}

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

static double metaDouble(const std::string &path, const char *key, double dflt)
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

/// Same padding scheme as CbcRootFixtureDump.hpp / mip-root-replay.cpp: drop
/// the redundant final row a padded capture carries.
static bool dropPadRow(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string meta = stem + ".meta";
  if (metaInt(meta, "paddedColumns", 0) <= 0)
    return false;
  const long capturedRows = metaInt(meta, "rows", -1);
  if (capturedRows < 0 || si.getNumRows() != (int)capturedRows + 1) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: meta says padded but rows=%d against captured %ld\n",
        baseName(stem).c_str(), si.getNumRows(), capturedRows);
    return false;
  }
  const int last = si.getNumRows() - 1;
  si.deleteRows(1, &last);
  return true;
}

static int restoreColTypes(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string path = stem + ".ctype";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: no .ctype sidecar\n", baseName(stem).c_str());
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
      fclose(fp);
      return -1;
    }
    if (si.isContinuous(idx)) {
      si.setInteger(idx);
      ++restored;
    }
  }
  fclose(fp);
  si.getColType(true);
  return restored;
}

/// Load the dense incumbent sidecar into `values` (resized to numberColumns).
static bool loadIncumbent(const std::string &stem, int numberColumns,
  std::vector< double > &values, double &incumbentObj)
{
  const std::string path = stem + ".incumbent";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    fprintf(stderr, "ERROR: %s: no .incumbent sidecar\n", baseName(stem).c_str());
    return false;
  }
  int sidecarCols = -1;
  if (fscanf(fp, "cols %d\n", &sidecarCols) != 1 || sidecarCols != numberColumns) {
    fprintf(stderr, "ERROR: %s: .incumbent is for %d columns, model has %d\n",
      baseName(stem).c_str(), sidecarCols, numberColumns);
    fclose(fp);
    return false;
  }
  if (fscanf(fp, "obj %lf\n", &incumbentObj) != 1) {
    fclose(fp);
    return false;
  }
  values.assign(numberColumns, 0.0);
  int idx = 0;
  double v = 0.0;
  while (fscanf(fp, "%d %lf\n", &idx, &v) == 2) {
    if (idx < 0 || idx >= numberColumns) {
      fclose(fp);
      return false;
    }
    values[idx] = v;
  }
  fclose(fp);
  return true;
}

/// Count integer columns where the LP relaxation's current value agrees
/// (within 10x primal tolerance) with the incumbent -- the same test
/// CbcHeuristicRINS::solution() itself uses to decide `nFix`, computed here
/// independently so --sweep can report it without instrumenting the class.
static int countExactAgreement(const OsiSolverInterface &si, const std::vector< double > &incumbent,
  double primalTolerance, int shallowDepth)
{
  const int n = si.getNumCols();
  const double *x = si.getColSolution();
  const double *lo = si.getColLower();
  const double *up = si.getColUpper();
  int nFix = 0;
  for (int j = 0; j < n; ++j) {
    if (si.isContinuous(j))
      continue;
    double value = incumbent[j];
    if (value < lo[j])
      value = lo[j];
    else if (value > up[j])
      value = up[j];
    const double nearest = floor(value + 0.5);
    if (fabs(x[j] - nearest) >= 10.0 * primalTolerance)
      continue;
    // Same shallowDepth_ filter CbcHeuristicRINS::solution() applies (modes
    // 1-3 additionally require the fixed value to be at/away from the
    // column's original lower bound; mode 3 also excludes columns used in a
    // prior solution, which this standalone count has no visibility into and
    // so treats as "not used" -- a slight over-count for mode 3 only).
    bool fix = false;
    switch (shallowDepth) {
    case 0: fix = true; break;
    case 1: fix = (nearest == lo[j]); break;
    case 2: fix = (nearest != lo[j]); break;
    case 3: fix = (nearest == lo[j]); break;
    default: fix = true; break;
    }
    if (fix)
      ++nFix;
  }
  return nFix;
}

/// Run CglProbing on a solver whose bounds already reflect RINS's fixing
/// (i.e. after setColLower==setColUpper on the agreeing columns), and report
/// how many ADDITIONAL columns it manages to fix purely by bound propagation
/// -- this is the empirical answer to "would more probing let us fix more
/// variables before delegating to the sub-MIP solver". Non-destructive: works
/// on a clone.
static int probeAdditionalFixes(const OsiSolverInterface &fixedSolver)
{
  OsiSolverInterface *clone = fixedSolver.clone();
  const int n = clone->getNumCols();
  int before = 0;
  {
    const double *lo0 = clone->getColLower();
    const double *up0 = clone->getColUpper();
    for (int j = 0; j < n; ++j)
      if (lo0[j] == up0[j])
        ++before;
  }

  CglProbing probing;
  probing.setMaxPass(5);
  probing.setMaxProbe(1000);
  probing.setMaxLook(1000);
  probing.setRowCuts(0); // tightening only, no cutting planes needed for this report
  probing.setUsingObjective(false);
  clone->resolve();
  if (!clone->isProvenOptimal()) {
    delete clone;
    return 0; // fixing already made the sub-LP infeasible; not this report's concern
  }
  OsiCuts cuts;
  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  probing.generateCuts(*clone, cuts, info);

  // Bound tightening comes back as OsiColCut objects, not applied to the
  // solver directly (generateCuts takes a const reference) -- apply them
  // here to see the effect.
  const int nColCuts = cuts.sizeColCuts();
  for (int i = 0; i < nColCuts; ++i) {
    const OsiColCut *cc = cuts.colCutPtr(i);
    const CoinPackedVector &lbs = cc->lbs();
    const CoinPackedVector &ubs = cc->ubs();
    for (int k = 0; k < lbs.getNumElements(); ++k)
      clone->setColLower(lbs.getIndices()[k], lbs.getElements()[k]);
    for (int k = 0; k < ubs.getNumElements(); ++k)
      clone->setColUpper(ubs.getIndices()[k], ubs.getElements()[k]);
  }

  int after = 0;
  {
    const double *lo1 = clone->getColLower();
    const double *up1 = clone->getColUpper();
    for (int j = 0; j < n; ++j)
      if (lo1[j] == up1[j])
        ++after;
  }
  delete clone;
  return after - before;
}

struct RunResult {
  bool found = false;
  double objValue = 0.0;
  double improvement = 0.0; // incumbentObj - objValue (assuming minimization
                             // sense as reported); negative means worse
  double timeSec = 0.0;
  int nFixExact = 0;
  int numberIntegers = 0;
  int probeExtra = -1; // -1 when not requested
};

/// One sweep point: build a fresh CbcModel around the fixture, attach one
/// heuristic instance configured with the given knobs, call solution() once.
static RunResult runOnePoint(const std::string &mps, const std::string &stem,
  const std::string &method, int shallow, double fixClose, int nodes, double fracSmall,
  double secBudget, int logLevel, bool probeReport, bool quiet)
{
  RunResult result;

  OsiClpSolverInterface si;
  ClpSimplex *lp = si.getModelPtr();
  lp->setLogLevel(0);
  if (si.readMps(mps.c_str())) {
    fprintf(stderr, "ERROR: failed to read %s\n", mps.c_str());
    return result;
  }
  dropPadRow(si, stem, quiet);
  restoreColTypes(si, stem, quiet);

  const std::string bas = stem + ".bas";
  bool haveBasis = false;
  if (fileExists(bas)) {
    if (lp->readBasis(bas.c_str()) >= 0) {
      haveBasis = true;
      si.setWarmStart(NULL);
    }
  }
  lp->setPerturbation(50);
  if (haveBasis) {
    si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
    si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
    si.resolve();
  } else {
    si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);
    si.initialSolve();
  }
  if (!si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: LP not optimal after warm start (%s)\n", stem.c_str());
    return result;
  }

  std::vector< double > incumbent;
  double incumbentObj = 0.0;
  if (!loadIncumbent(stem, si.getNumCols(), incumbent, incumbentObj))
    return result;

  double primalTolerance = 1.0e-7;
  si.getDblParam(OsiPrimalTolerance, primalTolerance);
  result.nFixExact = countExactAgreement(si, incumbent, primalTolerance, shallow);

  CbcModel model(si);
  model.setLogLevel(0);
  model.messageHandler()->setLogLevel(logLevel);
  model.setMaximumSeconds(secBudget);
  // Normally set up by branchAndBound() before any heuristic runs; needed
  // here since this tool calls the heuristic directly without ever entering
  // the tree -- without it, model.object(i) is null and RINS/VND crash.
  model.findIntegers(true);
  result.numberIntegers = model.numberIntegers();

  // Install the incumbent as the model's best solution: setBestSolution()
  // does NOT bump numberSolutions_ itself (it only exists so callers can
  // decide), so bump it explicitly -- this is exactly the state
  // CbcHeuristicRINS::solution()'s `numberSolutions_ < getSolutionCount()`
  // gate checks on its first call.
  model.setBestSolution(incumbent.data(), (int)incumbent.size(), incumbentObj, false);
  model.setSolutionCount(1);

  CbcHeuristic *heuristic = NULL;
  CbcHeuristicRINS *rins = NULL;
  if (method == "rins") {
    rins = new CbcHeuristicRINS(model);
    rins->setShallowDepth(shallow);
    rins->setFixCloseMaxDist(fixClose);
    rins->setNumberNodes(nodes);
    rins->setFractionSmall(fracSmall);
    heuristic = rins;
  } else if (method == "vnd") {
    heuristic = new CbcHeuristicVND(model);
    heuristic->setShallowDepth(shallow);
    heuristic->setNumberNodes(nodes);
    heuristic->setFractionSmall(fracSmall);
  } else {
    fprintf(stderr, "ERROR: unknown --method=%s (expected rins or vnd)\n", method.c_str());
    return result;
  }
  heuristic->setModel(&model);

  if (probeReport && rins) {
    // Replicate RINS's own exact-agreement fixing pass on a scratch clone,
    // purely to measure how much MORE probing could fix -- does not affect
    // the actual heuristic call below.
    OsiSolverInterface *scratch = si.clone();
    const int numberIntegers = model.numberIntegers();
    const int *integerVariable = model.integerVariable();
    const double *x = si.getColSolution();
    for (int i = 0; i < numberIntegers; ++i) {
      const int iColumn = integerVariable[i];
      double value = incumbent[iColumn];
      const double lo = scratch->getColLower()[iColumn];
      const double up = scratch->getColUpper()[iColumn];
      if (value < lo)
        value = lo;
      else if (value > up)
        value = up;
      const double nearest = floor(value + 0.5);
      if (fabs(x[iColumn] - nearest) < 10.0 * primalTolerance) {
        scratch->setColLower(iColumn, nearest);
        scratch->setColUpper(iColumn, nearest);
      }
    }
    result.probeExtra = probeAdditionalFixes(*scratch);
    delete scratch;
  }

  std::vector< double > newSolution(si.getNumCols(), 0.0);
  double heuristicValue = model.getCutoff();

  const double t0 = wallClock();
  const int rc = heuristic->solution(heuristicValue, newSolution.data());
  result.timeSec = wallClock() - t0;

  if (rc > 0) {
    result.found = true;
    result.objValue = heuristicValue;
    result.improvement = incumbentObj - heuristicValue; // positive = better
  }

  delete heuristic;
  return result;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <instance-or-fixture> [options]\n"
    "  --fixture-dir=PATH   default $CBC_RINS_FIXTURE_DIR or "
    "~/instances/mip-sanity-data/rinsFixtures\n"
    "  --method=rins|vnd    default rins\n"
    "  --shallow=N          default 0\n"
    "  --fix-close=X        default 0.4 (RINS only)\n"
    "  --nodes=N            default 200\n"
    "  --frac-small=X       default 0.5\n"
    "  --sec=N              default 30\n"
    "  --probe-report       report additional fixes CglProbing would find\n"
    "  --sweep              run a built-in shallow x fix-close x nodes grid\n"
    "  --log=N              default 0\n"
    "  --quiet\n",
    prog);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    usage(argv[0]);
    return 2;
  }

  std::string arg1;
  std::string fixtureDir;
  {
    const char *env = getenv("CBC_RINS_FIXTURE_DIR");
    const char *home = getenv("HOME");
    fixtureDir = (env && *env) ? env
                                : std::string(home ? home : ".") + "/instances/mip-sanity-data/rinsFixtures";
  }
  std::string method = "rins";
  int shallow = 0, nodes = 200, logLevel = 0;
  double fixClose = 0.4, fracSmall = 0.5, sec = 30.0;
  bool sweep = false, probeReport = false, quiet = false;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--fixture-dir=", 0) == 0)
      fixtureDir = a.substr(14);
    else if (a.rfind("--method=", 0) == 0)
      method = a.substr(9);
    else if (a.rfind("--shallow=", 0) == 0)
      shallow = atoi(a.c_str() + 10);
    else if (a.rfind("--fix-close=", 0) == 0)
      fixClose = atof(a.c_str() + 12);
    else if (a.rfind("--nodes=", 0) == 0)
      nodes = atoi(a.c_str() + 8);
    else if (a.rfind("--frac-small=", 0) == 0)
      fracSmall = atof(a.c_str() + 13);
    else if (a.rfind("--sec=", 0) == 0)
      sec = atof(a.c_str() + 6);
    else if (a == "--probe-report")
      probeReport = true;
    else if (a == "--sweep")
      sweep = true;
    else if (a.rfind("--log=", 0) == 0)
      logLevel = atoi(a.c_str() + 6);
    else if (a == "--quiet")
      quiet = true;
    else if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a.empty() || a[0] == '-') {
      fprintf(stderr, "Unknown option: %s\n", a.c_str());
      usage(argv[0]);
      return 2;
    } else if (arg1.empty())
      arg1 = a;
  }
  if (arg1.empty()) {
    usage(argv[0]);
    return 2;
  }

  const std::string stem = fixtureStem(arg1, fixtureDir);
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string name = baseName(stem);
  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no fixture problem file for stem %s (looked for %s)\n",
      stem.c_str(), mps.c_str());
    return 1;
  }

  if (shallow != 0 && !quiet) {
    fprintf(stderr,
      "NOTE: --shallow=%d requested; a fresh CbcModel treats THIS node's "
      "(branching-tightened) bounds as 'original', not the root's -- see "
      "CbcHeuristicRinsFixtureDump.hpp's caveat. shallow=0 is unaffected.\n",
      shallow);
  }

  if (!sweep) {
    const RunResult r = runOnePoint(mps, stem, method, shallow, fixClose, nodes, fracSmall,
      sec, logLevel, probeReport, quiet);
    printf("[rins-bench] %s: method=%s shallow=%d fixClose=%g nodes=%d fracSmall=%g "
           "nFixExact=%d/%d found=%d obj=%s improve=%s time=%.4fs",
      name.c_str(), method.c_str(), shallow, fixClose, nodes, fracSmall, r.nFixExact,
      r.numberIntegers, (int)r.found, r.found ? std::to_string(r.objValue).c_str() : "n/a",
      r.found ? std::to_string(r.improvement).c_str() : "n/a", r.timeSec);
    if (probeReport)
      printf(" probeExtraFixes=%d", r.probeExtra);
    printf("\n");
    return 0;
  }

  // Built-in sweep grid.
  static const int shallowGrid[] = { 0, 1, 2, 3 };
  static const double fixCloseGrid[] = { 0.0, 0.2, 0.4, 0.6 };
  static const int nodesGrid[] = { 50, 200, 1000 };

  if (!quiet) {
    fprintf(stderr,
      "NOTE: --sweep runs all 48 grid points back-to-back inside this one process.\n"
      "Confirmed empirically (2026-09) that some grid points -- especially later ones,\n"
      "and VND in particular -- can silently return a worse/no-solution result here than\n"
      "an isolated single-point run with the IDENTICAL options would (root cause not yet\n"
      "found; some state is leaking between successive CbcModel/OsiSolverInterface builds\n"
      "in one process). Treat --sweep as a fast, DIRECTIONAL ranking tool only -- always\n"
      "re-confirm any promising cell with a separate one-shot invocation (no --sweep)\n"
      "before trusting its exact found/improve value.\n");
  }
  printf("%-30s %-6s %-4s %-8s %-6s %-9s %-9s %-6s %-14s %-9s\n", "instance", "method",
    "shl", "fixClose", "nodes", "nFixExact", "numInt", "found", "improve", "time(s)");
  for (size_t si2 = 0; si2 < sizeof(shallowGrid) / sizeof(shallowGrid[0]); ++si2) {
    for (size_t ci = 0; ci < sizeof(fixCloseGrid) / sizeof(fixCloseGrid[0]); ++ci) {
      for (size_t ni = 0; ni < sizeof(nodesGrid) / sizeof(nodesGrid[0]); ++ni) {
        const RunResult r = runOnePoint(mps, stem, method, shallowGrid[si2], fixCloseGrid[ci],
          nodesGrid[ni], fracSmall, sec, logLevel, false, /*quiet=*/true);
        printf("%-30s %-6s %-4d %-8g %-6d %-9d %-9d %-6d %-14s %-9.4f\n", name.c_str(),
          method.c_str(), shallowGrid[si2], fixCloseGrid[ci], nodesGrid[ni], r.nFixExact,
          r.numberIntegers, (int)r.found,
          r.found ? std::to_string(r.improvement).c_str() : "n/a", r.timeSec);
      }
    }
  }

  return 0;
}
