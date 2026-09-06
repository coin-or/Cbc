/**
 * mip-root-replay — replay CBC's root-node processing (cut generation +
 * heuristics, plus a handful of initial B&B nodes) from a fixture, without
 * paying for pre-processing or the root LP solve again.
 *
 * @file mip-root-replay.cpp
 * @brief fast iteration loop for root cut-generation / heuristics strategies
 *
 * Loads a fixture written by CbcRootFixtureDump.hpp (CBC_DUMP_ROOT_FIXTURE) --
 * the preprocessed problem plus its optimal root LP basis -- warm-starts the
 * LP (0 iterations on a faithful fixture), builds a fresh CbcModel directly
 * around it, attaches the same default cut generators / heuristics
 * CbcStrategyDefault would give the normal `cbc` command line, and calls
 * branchAndBound() with a small node limit. Since the fixture already sits at
 * the exact point CBC would reach after pre-processing + the root LP, this
 * reaches root cut generation / heuristics in milliseconds instead of a full
 * solve -- see ROOT-FIXTURES.md for the fixture format and generation
 * workflow (`gen-root-fixtures`).
 *
 * Deliberately does NOT re-run pre-processing: CbcStrategyDefault is used
 * without calling setupPreProcessing(), so strategy_->setupOther() only sets
 * numberStrong_/numberBeforeTrust_ and does not touch the solver. The
 * fixture's solver is already what CBC's own preprocessed model looked like.
 *
 * Usage:
 *   mip-root-replay <instance-name> [options]     (looks up the fixture in
 *                                                   --fixture-dir / mip-sanity-data)
 *   mip-root-replay <stem-or-any-fixture-file> [options]
 *
 * Options:
 *   --fixture-dir=PATH  Fixture directory (default: $CBC_ROOT_FIXTURE_DIR or
 *                        ~/instances/mip-sanity-data/rootFixtures)
 *   --tag=TAG            Fixture tag (default "root")
 *   --nodes=N            Max nodes to explore, root counts as node 1 (default 16)
 *   --sec=N               Wall time backstop for the replay itself (default 60)
 *   --no-cuts             Skip attaching cut generators (heuristics only)
 *   --no-heur             Skip attaching heuristics (cuts only)
 *   --log=N               CBC log level (default 0)
 *   --data-dir=PATH       mip-sanity-data checkout, for bks.tsv gap reporting
 *                         (default: ../mip-sanity-data relative to this binary)
 *   --quiet               Suppress warnings (e.g. missing sidecars)
 *
 * Exit codes: 0 loaded and ran; 1 fixture load failed; 2 usage error.
 */

#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "ClpSimplex.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

static double wallClock()
{
  return CoinGetTimeOfDay();
}

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

/// Reduce any fixture file name (or a bare instance name) to the shared stem.
static std::string fixtureStem(const std::string &arg, const std::string &fixtureDir,
  const std::string &tag)
{
  static const char *suffixes[] = { ".mps.gz", ".mps", ".bas", ".sol", ".ctype", ".meta",
    ".dumplog" };
  for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
    const std::string suf(suffixes[i]);
    if (arg.size() > suf.size() && arg.compare(arg.size() - suf.size(), suf.size(), suf) == 0)
      return arg.substr(0, arg.size() - suf.size());
  }
  // A bare name (no path separator, no known suffix): look it up in the
  // fixture directory as "<name>.<tag>".
  if (arg.find('/') == std::string::npos)
    return fixtureDir + "/" + arg + "." + tag;
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

/// Same padding scheme as CbcRootFixtureDump.hpp / CbcClqFixtureDump.hpp: drop
/// the redundant final row a padded capture carries, so the model matches
/// what was captured.
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

/// See CbcClqFixtureDump.hpp / bkclique-bench.cpp's restoreColTypes for full
/// rationale: MPS cannot express "fixed and integer", so a `.ctype` sidecar
/// carries integrality for columns pre-processing fixed.
static int restoreColTypes(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string path = stem + ".ctype";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: no .ctype sidecar; integer columns fixed at "
                      "capture will read back continuous\n",
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
  si.getColType(true); // recompute the cached Binary/GeneralInteger derivation
  return restored;
}

static std::vector< std::string > splitTab(const std::string &line)
{
  std::vector< std::string > fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t'))
    fields.push_back(field);
  return fields;
}

static std::vector< std::string > lookupRow(const std::string &tsvPath, const std::string &key)
{
  std::ifstream in(tsvPath);
  if (!in.is_open())
    return {};
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (first) {
      first = false;
      continue;
    }
    if (line.empty())
      continue;
    std::vector< std::string > fields = splitTab(line);
    if (!fields.empty() && fields[0] == key)
      return fields;
  }
  return {};
}

/// Thin wrapper so --no-cuts/--no-heur can skip either half of
/// CbcStrategyDefault's setup without reimplementing it.
class ReplayStrategy : public CbcStrategyDefault {
public:
  ReplayStrategy(bool doCuts, bool doHeur, int numberStrong, int numberBeforeTrust)
    : CbcStrategyDefault(1, numberStrong, numberBeforeTrust)
    , doCuts_(doCuts)
    , doHeur_(doHeur)
  {
  }
  virtual CbcStrategy *clone() const { return new ReplayStrategy(*this); }
  virtual void setupCutGenerators(CbcModel &model)
  {
    if (doCuts_)
      CbcStrategyDefault::setupCutGenerators(model);
  }
  virtual void setupHeuristics(CbcModel &model)
  {
    if (doHeur_)
      CbcStrategyDefault::setupHeuristics(model);
  }

private:
  bool doCuts_, doHeur_;
};

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <instance-name-or-fixture-stem> [options]\n"
    "  --fixture-dir=PATH  fixture directory\n"
    "                      (default: $CBC_ROOT_FIXTURE_DIR or "
    "~/instances/mip-sanity-data/rootFixtures)\n"
    "  --tag=TAG           fixture tag (default \"root\")\n"
    "  --nodes=N            max nodes, root counts as node 1 (default 16)\n"
    "  --sec=N              wall time backstop for the replay (default 60)\n"
    "  --no-cuts            skip attaching cut generators\n"
    "  --no-heur            skip attaching heuristics\n"
    "  --log=N              CBC log level (default 0)\n"
    "  --data-dir=PATH      mip-sanity-data checkout, for bks.tsv gap reporting\n"
    "  --quiet              suppress warnings\n",
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
    const char *env = getenv("CBC_ROOT_FIXTURE_DIR");
    const char *home = getenv("HOME");
    fixtureDir = (env && *env) ? env
                                : std::string(home ? home : ".") + "/instances/mip-sanity-data/rootFixtures";
  }
  std::string tag = "root";
  int nodes = 16;
  double sec = 60.0;
  bool doCuts = true, doHeur = true, quiet = false;
  int logLevel = 0;
  std::string dataDir;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("--fixture-dir=", 0) == 0)
      fixtureDir = a.substr(14);
    else if (a.rfind("--tag=", 0) == 0)
      tag = a.substr(6);
    else if (a.rfind("--nodes=", 0) == 0)
      nodes = atoi(a.c_str() + 8);
    else if (a.rfind("--sec=", 0) == 0)
      sec = atof(a.c_str() + 6);
    else if (a == "--no-cuts")
      doCuts = false;
    else if (a == "--no-heur")
      doHeur = false;
    else if (a.rfind("--log=", 0) == 0)
      logLevel = atoi(a.c_str() + 6);
    else if (a.rfind("--data-dir=", 0) == 0)
      dataDir = a.substr(11);
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

  const std::string stem = fixtureStem(arg1, fixtureDir, tag);
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";
  const std::string name = baseName(stem);

  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no fixture problem file for stem %s (looked for %s)\n",
      stem.c_str(), mps.c_str());
    return 1;
  }

  OsiClpSolverInterface si;
  ClpSimplex *lp = si.getModelPtr();
  lp->setLogLevel(0);

  if (si.readMps(mps.c_str())) {
    fprintf(stderr, "ERROR: failed to read %s\n", mps.c_str());
    return 1;
  }

  dropPadRow(si, stem, quiet);
  restoreColTypes(si, stem, quiet);

  bool haveBasis = false;
  if (fileExists(bas)) {
    if (lp->readBasis(bas.c_str()) < 0) {
      fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
    } else {
      haveBasis = true;
      // Push the freshly-read status into OsiClp's cached basis_, or resolve()
      // below installs the stale cache over it. See bkclique-bench.cpp's
      // loadFixture for the full rationale.
      si.setWarmStart(NULL);
    }
  } else if (!quiet) {
    fprintf(stderr, "WARNING: no basis %s; solving cold\n", bas.c_str());
  }

  lp->setPerturbation(50);

  const double t0 = wallClock();
  if (haveBasis) {
    si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
    si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
    si.resolve();
  } else {
    si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);
    si.initialSolve();
  }
  const double loadTime = wallClock() - t0;
  const int warmStartIters = si.getIterationCount();

  if (!si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: LP not optimal after warm start (%s)\n", stem.c_str());
    return 1;
  }
  if (haveBasis && warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured "
                    "basis did not survive, so this replay starts from a "
                    "different vertex\n",
      name.c_str(), warmStartIters);
  }

  printf("[replay] %s: rows=%d cols=%d warmStartIters=%d loadTime=%.4fs\n",
    name.c_str(), si.getNumRows(), si.getNumCols(), warmStartIters, loadTime);

  CbcModel model(si);
  model.setLogLevel(0);
  model.messageHandler()->setLogLevel(logLevel);
  model.setMaximumNodes(nodes);
  model.setMaximumSeconds(sec);

  // Same recipe the normal `cbc` command line uses for its default cut
  // generators / heuristics (see CbcSolver.cpp's babExecuteSearchAndPostprocess,
  // "CbcStrategyDefault strategy(1, babModel_->numberStrong(),
  // babModel_->numberBeforeTrust())"). Deliberately not calling
  // setupPreProcessing(): the fixture's solver is already preprocessed, and
  // CbcStrategyDefault::setupOther() only preprocesses when that is requested,
  // so it is a no-op here beyond setting numberStrong_/numberBeforeTrust_.
  ReplayStrategy strategy(doCuts, doHeur, model.numberStrong(), model.numberBeforeTrust());
  model.setStrategy(strategy);

  const double t1 = wallClock();
  model.branchAndBound();
  const double bbTime = wallClock() - t1;

  const bool optimal = model.isProvenOptimal();
  const bool infeasible = model.isProvenInfeasible();
  const double bound = model.getBestPossibleObjValue();
  const bool haveSolution = model.getSolutionCount() > 0 || (optimal && !infeasible);
  const double best = haveSolution ? model.getObjValue() : COIN_DBL_MAX;
  double gapPct = -1.0;
  if (haveSolution && std::fabs(best) > 1e-10)
    gapPct = 100.0 * std::fabs(best - bound) / std::fabs(best);

  printf("[replay] %s: nodes=%lld optimal=%d infeasible=%d bound=%.10g best=%s gap=%s "
         "bbTime=%.4fs\n",
    name.c_str(), (long long)model.getNodeCount(), (int)optimal, (int)infeasible, bound,
    haveSolution ? std::to_string(best).c_str() : "none",
    gapPct >= 0 ? (std::to_string(gapPct) + "%").c_str() : "n/a", bbTime);

  if (dataDir.empty()) {
    // Default: ../mip-sanity-data relative to this binary's own build
    // directory (test/), so it works out of the box from an in-tree build.
    dataDir = "mip-sanity-data";
  }
  const std::vector< std::string > bksRow = lookupRow(dataDir + "/bks.tsv", name);
  if (!bksRow.empty() && bksRow.size() > 1) {
    const double bks = atof(bksRow[1].c_str());
    printf("[replay] %s: bks=%.10g\n", name.c_str(), bks);
  }

  return 0;
}
