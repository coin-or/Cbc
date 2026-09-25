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
#include "CbcParameters.hpp"
#include "CbcSolverCutSetup.hpp"
#include "CbcSolverHeuristics.hpp"
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

/**
 * Load a `.debugsol` sidecar written by cbcRootFixtureWriteDebugSol() and
 * attach it to `si` via activateRowCutDebugger(), so every debugger->
 * invalidCut()/CoinAssert check already wired throughout CbcModel.cpp /
 * CbcCutGenerator.cpp fires exactly as it would on the real solve -- turning
 * an invalid-cut bug (a cut generator silently cutting off the true optimum)
 * into an instant, offline-reproducible failure instead of a full solve.
 *
 * Returns true if a debugger was attached. Absence of the sidecar is silent,
 * not a warning: it only exists when the fixture was captured with
 * "-debugCuts <file>" active, which is the exception, not the rule.
 */
static bool loadDebugSol(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string path = stem + ".debugsol";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return false;

  int sidecarCols = -1;
  if (fscanf(fp, "cols %d\n", &sidecarCols) != 1 || sidecarCols != si.getNumCols()) {
    fprintf(stderr, "ERROR: %s: .debugsol is for %d columns, model has %d; ignoring it\n",
      baseName(stem).c_str(), sidecarCols, si.getNumCols());
    fclose(fp);
    return false;
  }

  std::vector< double > values(sidecarCols, 0.0);
  int idx = 0;
  double value = 0.0;
  int loaded = 0, nonzero = 0;
  while (fscanf(fp, "%d %lf\n", &idx, &value) == 2) {
    if (idx < 0 || idx >= sidecarCols) {
      fprintf(stderr, "ERROR: %s: .debugsol names column %d, out of range\n",
        baseName(stem).c_str(), idx);
      fclose(fp);
      return false;
    }
    values[idx] = value;
    ++loaded;
    if (value != 0.0)
      ++nonzero;
  }
  fclose(fp);

  // activateRowCutDebugger(solution) builds a fresh OsiRowCutDebugger from
  // this array against si's CURRENT columns -- exactly the state the sidecar
  // was written in, since it was captured in the same preprocessed column
  // order as the .mps.gz this si was just loaded from.
  si.activateRowCutDebugger(values.data());
  if (!quiet)
    fprintf(stderr, "[replay] %s: loaded debug/reference solution (%d values, "
                    "%d nonzero) -- cuts will be checked against it\n",
      baseName(stem).c_str(), loaded, nonzero);
  return true;
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

static std::vector< std::string > lookupRow(const std::string &tsvPath, const std::string &key,
  std::vector< std::string > *headerOut = nullptr)
{
  std::ifstream in(tsvPath);
  if (!in.is_open())
    return {};
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (first) {
      first = false;
      if (headerOut)
        *headerOut = splitTab(line);
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

/// Thin wrapper so --no-cuts/--no-heur can skip either half of root setup.
///
/// Cut generators are installed via the *real* installCutGenerators()
/// (CbcSolverCutSetup.hpp) -- the exact free function the production `cbc`
/// CLI calls from CbcSolver::configureCutGenerators() -- rather than
/// CbcStrategyDefault::setupCutGenerators(). The latter looked like the
/// obvious choice (it's what a hand-built CbcStrategy would normally use),
/// but CbcStrategy.cpp's setupCutGenerators() hardcodes `genFlags=0` in its
/// non-CBC_GENERATE_TEST build path, which silently adds *zero* cut
/// generators (dead code, apparently a leftover from CBC_GENERATE_TEST
/// scaffolding) -- confirmed by a replay that produced exactly one cut pass
/// with 0 tight cuts on an instance the real CLI cuts substantially. Calling
/// the real free function instead, with the same scalar defaults
/// CbcSolver::run() uses (complicatedInteger=0, dominatedCuts=false,
/// cgraphMode="on", oldCliqueMode=parameters' CLIQUECUTS default,
/// maxCallsBK=1000, bkClqExtMethod=4, bkPivotingStrategy=Weight,
/// oddWExtMethod=2, mixedRoundStrategy=1; see CbcSolver.cpp's `run()`
/// whereFrom==2 block) reproduces the exact generator set/tuning the CLI
/// installs, then reapplies the CLI's own minimumDrop /
/// maximumCutPassesAtRoot / maximumCutPasses formulas (CbcSolver.cpp ~8194-
/// 8211) so an unmodified replay matches the real CLI's root behavior
/// pass-for-pass. --pass-cuts=N (passCutsOverride) is applied last, so a
/// sweep can substitute a candidate tiering rule for the CLI's default one.
///
/// Heuristics are attached via the real doHeuristics() used by the `cbc`
/// command line itself (rounding, Feasibility Pump, RINS, diving,
/// Feasibility Jump, the FPump->FJ fallback, ...) rather than
/// CbcStrategyDefault::setupHeuristics()'s much smaller fixed set (rounding
/// only) -- this is what lets replay experiments faithfully compare against
/// full-CLI heuristic behavior.
/// Parse --rins-schedule=/--vnd-schedule= string values into
/// HeuristicScheduleMode. Empty string means "not requested" (legacy,
/// unchanged). Exits with an error message on an unrecognized value.
static HeuristicScheduleMode parseScheduleModeArg(const std::string &s)
{
  if (s.empty() || s == "legacy")
    return HeuristicScheduleMode::Legacy;
  if (s == "depth")
    return HeuristicScheduleMode::EveryKDepth;
  if (s == "nodes")
    return HeuristicScheduleMode::EveryKNodes;
  if (s == "nodes-no-improve")
    return HeuristicScheduleMode::EveryKNodesNoImprove;
  fprintf(stderr, "Unknown schedule mode: %s\n", s.c_str());
  exit(2);
}

class ReplayStrategy : public CbcStrategyDefault {
public:
  ReplayStrategy(bool doCuts, bool doHeur, int numberStrong, int numberBeforeTrust,
    CbcParameters &cutParams, CbcParameters *heurParams = NULL, int passCutsOverride = 0,
    double minDropScale = 1.0, const std::string &rinsSchedule = "", int rinsScheduleK = 1,
    const std::string &vndSchedule = "", int vndScheduleK = 1)
    : CbcStrategyDefault(1, numberStrong, numberBeforeTrust)
    , doCuts_(doCuts)
    , doHeur_(doHeur)
    , cutParams_(cutParams)
    , heurParams_(heurParams)
    , passCutsOverride_(passCutsOverride)
    , minDropScale_(minDropScale)
    , rinsSchedule_(rinsSchedule)
    , rinsScheduleK_(rinsScheduleK)
    , vndSchedule_(vndSchedule)
    , vndScheduleK_(vndScheduleK)
  {
  }
  virtual CbcStrategy *clone() const { return new ReplayStrategy(*this); }
  virtual void setupCutGenerators(CbcModel &model)
  {
    if (!doCuts_)
      return;
    installCutGenerators(model, cutParams_, /*complicatedInteger=*/0,
      /*dominatedCuts=*/false, /*cgraphMode=*/"on",
      /*oldCliqueMode=*/cutParams_[CbcParam::CLIQUECUTS]->modeVal(),
      /*maxCallsBK=*/1000, /*bkClqExtMethod=*/4,
      CoinBronKerbosch::PivotingStrategy::Weight,
      /*oddWExtMethod=*/2, /*mixedRoundStrategy=*/1);
    // Same recipe as CbcSolver.cpp's run() (whereFrom==2 block): minimum drop
    // scaled off the root objective, then the CLI's cols-tiered rule with an
    // OR-widened top tier: numCols<500 OR numRows<500 -> minDrop ignored, up
    // to 100 passes; numCols<5000 -> minDrop-limited, up to 100 passes; else
    // -> minDrop-limited, up to 50 passes. Overridden by --pass-cuts=N /
    // --min-drop-scale=X when passed explicitly (minDropScale_ multiplies
    // minimumDrop even in the default path, so --min-drop-scale=X alone can
    // still be used to explore alternative minDrop thresholds without a
    // --pass-cuts override).
    double minimumDrop = fabs(model.solver()->getObjValue()) * 1.0e-5 + 1.0e-5;
    minimumDrop = std::min(5.0e-2, minimumDrop);
    if (passCutsOverride_ != 0) {
      model.setMaximumCutPassesAtRoot(passCutsOverride_);
    } else {
      int numCols = model.getNumCols();
      int numRows = model.getNumRows();
      if (numCols < 500 || numRows < 500)
        model.setMaximumCutPassesAtRoot(-100);
      else if (numCols < 5000)
        model.setMaximumCutPassesAtRoot(100);
      else
        model.setMaximumCutPassesAtRoot(50);
    }
    model.setMinimumDrop(minimumDrop * minDropScale_);
    model.setMaximumCutPasses(4);
  }
  virtual void setupHeuristics(CbcModel &model)
  {
    if (!doHeur_)
      return;
    if (heurParams_)
      doHeuristics(&model, 1, *heurParams_, /*noPrinting_=*/1, /*initialPumpTune=*/0);
    else
      CbcStrategyDefault::setupHeuristics(model);
    // Apply --rins-schedule=/--vnd-schedule= overrides here, right after the
    // heuristics doHeuristics()/CbcStrategyDefault attached are actually
    // present on the model -- this runs synchronously inside
    // branchAndBound(), before the tree search starts, so it is the correct
    // (and only) place to reach the RINS/VND instances before they fire.
    if (rinsSchedule_.empty() && vndSchedule_.empty())
      return;
    const HeuristicScheduleMode rinsMode = parseScheduleModeArg(rinsSchedule_);
    const HeuristicScheduleMode vndMode = parseScheduleModeArg(vndSchedule_);
    for (int i = 0; i < model.numberHeuristics(); ++i) {
      CbcHeuristic *h = model.heuristic(i);
      if (!h)
        continue;
      const std::string hn = h->heuristicName() ? h->heuristicName() : "";
      if (!rinsSchedule_.empty() && hn == "RINS") {
        h->setScheduleMode(rinsMode, rinsScheduleK_);
      } else if (!vndSchedule_.empty() && hn == "VND") {
        h->setScheduleMode(vndMode, vndScheduleK_);
      }
    }
  }

private:
  bool doCuts_, doHeur_;
  CbcParameters &cutParams_;
  CbcParameters *heurParams_;
  int passCutsOverride_;
  double minDropScale_;
  std::string rinsSchedule_;
  int rinsScheduleK_;
  std::string vndSchedule_;
  int vndScheduleK_;
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
    "  --quiet              suppress warnings\n"
    "  --minimal-heur       use CbcStrategyDefault's bare rounding-only "
    "heuristic set\n"
    "                       instead of the real cbc CLI default set "
    "(doHeuristics())\n"
    "  --fpump=on|off       toggle Feasibility Pump (default: on, the real "
    "CLI default)\n"
    "  --fj=off|on|before|both   toggle Feasibility Jump (default: on, the "
    "real CLI default)\n"
    "  --fj-after-fpump=0|1|2 FJ fallback seeded from FPump's failed attempt; "
    "2=fallback-only,\n"
    "                       FJ never runs standalone (default 2)\n"
    "  --fj-effort=N        fixed FJ effort budget (default 0 = NNZ-scaled)\n"
    "  --fj-effort-mult=N   NNZ multiplier for FJ effort budget (default 1024)\n"
    "  --fj-stall=N         NNZ multiplier for FJ stall termination (default 256)\n"
    "  --fj-max-sol=N       stop FJ after this many solutions per call (default 1)\n"
    "  --fj-only-no-sol=0|1 only run FJ while no incumbent exists (default 1)\n"
    "  --fj-max-calls=N     cap on total FJ invocations (default 0 = unlimited)\n"
    "  --fj-depth=N         run FJ every N tree levels (default 0 = root only)\n"
    "  --jump-root-places=STR  jumpRootPlaces override, e.g. \"L\" (letters "
    "from LcCP; default: unset, uses CbcHeuristic's built-in \"LcC\")\n"
    "  --pump-root-places=STR  pumpRootPlaces override, e.g. \"L\" (default: "
    "unset, uses the legacy pumpTune-driven after-cuts logic)\n"
    "  --pass-cuts=N        override maximumCutPassesAtRoot after strategy setup\n"
    "                       (default: leave CbcStrategyDefault's own rule alone).\n"
    "                       Same encoding as the real CLI's -passCuts: positive N\n"
    "                       stops early once the minimum-drop test fails; negative\n"
    "                       N (abs value used as the pass cap) ignores minimum drop.\n"
    "  --min-drop-scale=X   multiply minimumDrop by X before applying it (default 1.0)\n"
    "  --gmi=off|root|ifmove|on       GMI cut mode (default: root, matching the\n"
    "                       real CLI's shipped default -- NOT CbcParameters()'s\n"
    "                       own out-of-the-box default of off)\n"
    "  --landp=off|root|ifmove|on     LandP cut mode (default: root, same rationale)\n"
    "  --redsplit2=off|root|ifmove|on RedSplit2 cut mode (default: root, same rationale)\n"
    "  --seed=N             CbcModel random seed (default: 1, matches CbcModel's\n"
    "                       own default). Vary this across repeats of the same\n"
    "                       config to average out branching tie-break noise\n"
    "                       (single-threaded B&B is otherwise deterministic).\n"
    "  --rins-schedule=legacy|depth|nodes|nodes-no-improve\n"
    "                       override RINS's re-invocation policy (default: legacy,\n"
    "                       the existing howOften_/decayFactor_ node-count math).\n"
    "                       See HeuristicScheduleMode in CbcHeuristic.hpp.\n"
    "  --rins-schedule-k=N  the K parameter for --rins-schedule (default 1)\n"
    "  --vnd-schedule=legacy|depth|nodes|nodes-no-improve   same, for VND\n"
    "  --vnd-schedule-k=N   the K parameter for --vnd-schedule (default 1)\n"
    "\n"
    "Invalid-cut / debug-cuts reproduction:\n"
    "  If <stem>.debugsol exists (written by CbcRootFixtureDump.hpp when the\n"
    "  dumping run was started with \"-debugCuts <solfile>\"), it is loaded\n"
    "  automatically and attached to the solver via activateRowCutDebugger()\n"
    "  before branchAndBound() runs. Every debugger->invalidCut()/CoinAssert\n"
    "  check already wired throughout CbcModel.cpp/CbcCutGenerator.cpp then\n"
    "  fires exactly as it would on the real solve, reproducing a false\n"
    "  \"proven infeasible\"/wrong-cut bug in milliseconds instead of a full\n"
    "  run. No flag is needed to opt in; there is no flag to opt out.\n",
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
  bool doCuts = true, doHeur = true, quiet = false, minimalHeur = false;
  int logLevel = 0;
  std::string dataDir;

  // Feasibility Jump / Feasibility Pump overrides, applied on top of
  // CbcParameters' own real CLI defaults (see below). NULL/unset (empty
  // string) means "leave the default alone".
  std::string fpumpMode; // "on"/"off"
  std::string fjMode; // "off"/"on"/"before"/"both"
  int fjAfterFPump = -1, fjEffort = -1, fjEffortMult = -1, fjStall = -1;
  int fjMaxSol = -1, fjOnlyNoSol = -1, fjMaxCalls = -1, fjDepth = -1;
  std::string jumpRootPlacesOverride, pumpRootPlacesOverride;
  int passCutsOverride = 0;
  double minDropScale = 1.0;
  // GMI/LandP/RedSplit2 cut-generator mode overrides. Default to "root" --
  // matching CbcSolver.cpp's initialize() shipped defaults (REDSPLIT2CUTS/
  // GMICUTS/LANDPCUTS all "root", see eeed3493) -- NOT the bare
  // CbcParameters() constructor default ("ifmove"/"off"/"off"), so an
  // unmodified replay (no --gmi=/--landp=/--redsplit2= flag) matches what
  // the real `cbc` CLI actually ships today, not CbcParameters' own
  // out-of-the-box defaults which predate that change.
  std::string gmiMode = "root", landpMode = "root", redsplit2Mode = "root";
  unsigned int randomSeed = 1; // matches CbcModel's own default
  std::string rinsScheduleMode, vndScheduleMode; // empty = legacy (no change)
  int rinsScheduleK = 1, vndScheduleK = 1;

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
    else if (a == "--minimal-heur")
      minimalHeur = true;
    else if (a.rfind("--fpump=", 0) == 0)
      fpumpMode = a.substr(8);
    else if (a.rfind("--fj=", 0) == 0)
      fjMode = a.substr(5);
    else if (a.rfind("--fj-after-fpump=", 0) == 0)
      fjAfterFPump = atoi(a.c_str() + 17);
    else if (a.rfind("--fj-effort-mult=", 0) == 0)
      fjEffortMult = atoi(a.c_str() + 17);
    else if (a.rfind("--fj-effort=", 0) == 0)
      fjEffort = atoi(a.c_str() + 12);
    else if (a.rfind("--fj-stall=", 0) == 0)
      fjStall = atoi(a.c_str() + 11);
    else if (a.rfind("--fj-max-sol=", 0) == 0)
      fjMaxSol = atoi(a.c_str() + 13);
    else if (a.rfind("--fj-only-no-sol=", 0) == 0)
      fjOnlyNoSol = atoi(a.c_str() + 17);
    else if (a.rfind("--fj-max-calls=", 0) == 0)
      fjMaxCalls = atoi(a.c_str() + 15);
    else if (a.rfind("--fj-depth=", 0) == 0)
      fjDepth = atoi(a.c_str() + 11);
    else if (a.rfind("--jump-root-places=", 0) == 0)
      jumpRootPlacesOverride = a.substr(19);
    else if (a.rfind("--pump-root-places=", 0) == 0)
      pumpRootPlacesOverride = a.substr(19);
    else if (a.rfind("--pass-cuts=", 0) == 0)
      passCutsOverride = atoi(a.c_str() + 12);
    else if (a.rfind("--min-drop-scale=", 0) == 0)
      minDropScale = atof(a.c_str() + 17);
    else if (a.rfind("--gmi=", 0) == 0)
      gmiMode = a.substr(6);
    else if (a.rfind("--landp=", 0) == 0)
      landpMode = a.substr(8);
    else if (a.rfind("--redsplit2=", 0) == 0)
      redsplit2Mode = a.substr(12);
    else if (a.rfind("--seed=", 0) == 0)
      randomSeed = (unsigned int)strtoul(a.c_str() + 7, NULL, 10);
    else if (a.rfind("--rins-schedule=", 0) == 0)
      rinsScheduleMode = a.substr(16);
    else if (a.rfind("--rins-schedule-k=", 0) == 0)
      rinsScheduleK = atoi(a.c_str() + 18);
    else if (a.rfind("--vnd-schedule=", 0) == 0)
      vndScheduleMode = a.substr(15);
    else if (a.rfind("--vnd-schedule-k=", 0) == 0)
      vndScheduleK = atoi(a.c_str() + 17);
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
  // The bare instance name, stripped of the ".<tag>" suffix fixtureStem()
  // appends (e.g. "pk1.root" -> "pk1") -- this is the key bks.tsv/limits.tsv
  // rows are actually indexed by, not the fixture stem's basename.
  const std::string instanceName = (name.size() > tag.size() + 1
                                      && name.compare(name.size() - tag.size() - 1, tag.size() + 1, "." + tag) == 0)
    ? name.substr(0, name.size() - tag.size() - 1)
    : name;

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
  // Must run before CbcModel model(si) below (whose copy constructor is what
  // carries rowCutDebugger_ into the model's own solver): ordering relative
  // to the LP resolve does not matter, activateRowCutDebugger() only sets up
  // the known-solution bookkeeping and does not touch the LP.
  const bool haveDebugSol = loadDebugSol(si, stem, quiet);

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
  // The LP relaxation bound *before any cuts* -- the fixture is captured at
  // exactly this point (preprocessed problem, optimal root basis, entering
  // the cut-generation loop for the first time) -- so this is the correct
  // baseline for a "gap closed by cuts" metric, as opposed to the raw
  // post-cuts bound alone which isn't comparable across instances.
  const double lpBound = si.getObjValue();

  printf("[replay] %s: rows=%d cols=%d warmStartIters=%d loadTime=%.4fs lpBound=%.10g "
         "debugSol=%d\n",
    name.c_str(), si.getNumRows(), si.getNumCols(), warmStartIters, loadTime, lpBound,
    (int)haveDebugSol);

  CbcModel model(si);
  model.setLogLevel(0);
  model.messageHandler()->setLogLevel(logLevel);
  model.setMaximumNodes(nodes);
  model.setMaximumSeconds(sec);
  model.setRandomSeed(randomSeed);

  // CbcParameters() self-initializes with the exact same defaults the real
  // `cbc` command line uses (CbcParameters::init() -> addCbcParams() ->
  // setDefaults()), independent of any CbcSolver instance -- so this is a
  // faithful baseline (Feasibility Pump on, Feasibility Jump off, etc.),
  // not a hand-picked subset. Overridden below only for the specific FJ/FPump
  // knobs this replay tool exposes.
  CbcParameters params;
  if (!fpumpMode.empty())
    params[CbcParam::FPUMP]->setKwdVal(fpumpMode);
  if (!fjMode.empty())
    params[CbcParam::FEASIBILITYJUMP]->setKwdVal(fjMode);
  if (fjAfterFPump >= 0)
    params[CbcParam::FEASIBILITYJUMPAFTERFPUMP]->setVal(fjAfterFPump);
  if (fjEffort >= 0)
    params[CbcParam::FEASIBILITYJUMPEFFORT]->setVal(fjEffort);
  if (fjEffortMult >= 0)
    params[CbcParam::FEASIBILITYJUMPEFFORTMULT]->setVal(fjEffortMult);
  if (fjStall >= 0)
    params[CbcParam::FEASIBILITYJUMPSTALL]->setVal(fjStall);
  if (fjMaxSol >= 0)
    params[CbcParam::FEASIBILITYJUMPMAXSOL]->setVal(fjMaxSol);
  if (fjOnlyNoSol >= 0)
    params[CbcParam::FEASIBILITYJUMPONLYNOSOL]->setVal(fjOnlyNoSol);
  if (fjMaxCalls >= 0)
    params[CbcParam::FEASIBILITYJUMPMAXCALLS]->setVal(fjMaxCalls);
  if (fjDepth >= 0)
    params[CbcParam::FEASIBILITYJUMPDEPTH]->setVal(fjDepth);
  if (!jumpRootPlacesOverride.empty())
    params[CbcParam::JUMPROOTPLACES]->setVal(jumpRootPlacesOverride);
  if (!pumpRootPlacesOverride.empty())
    params[CbcParam::PUMPROOTPLACES]->setVal(pumpRootPlacesOverride);
  // See default-init rationale above the flag parsing block: these three
  // default to "root" here (not CbcParameters()'s own "off"/"off"/"ifmove"
  // out-of-the-box defaults) to match CbcSolver.cpp's shipped CLI defaults.
  params[CbcParam::GMICUTS]->setVal(gmiMode);
  params[CbcParam::LANDPCUTS]->setVal(landpMode);
  params[CbcParam::REDSPLIT2CUTS]->setVal(redsplit2Mode);

  // Same recipe the normal `cbc` command line uses for its default cut
  // generators (see CbcSolver.cpp's babExecuteSearchAndPostprocess,
  // "CbcStrategyDefault strategy(1, babModel_->numberStrong(),
  // babModel_->numberBeforeTrust())"). Deliberately not calling
  // setupPreProcessing(): the fixture's solver is already preprocessed, and
  // CbcStrategyDefault::setupOther() only preprocesses when that is requested,
  // so it is a no-op here beyond setting numberStrong_/numberBeforeTrust_.
  // Heuristics are attached via doHeuristics() (the same function
  // CbcSolver::configureHeuristics() calls), for full parity with the real
  // CLI's default heuristic set, unless --minimal-heur asks for
  // CbcStrategyDefault's much smaller bare-rounding fallback instead.
  ReplayStrategy strategy(doCuts, doHeur, model.numberStrong(), model.numberBeforeTrust(),
    params, minimalHeur ? NULL : &params, passCutsOverride, minDropScale,
    rinsScheduleMode, rinsScheduleK, vndScheduleMode, vndScheduleK);
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
  // bks.tsv schema varies by collection: the simple MIPLIB-set format is
  // "instance\tbks", while mip-sanity-data's own bks.tsv is
  // "instance\tstatus\tobjective\tsense\tsource". Detect the right column
  // from the header instead of hardcoding index 1, so both work.
  std::vector< std::string > bksHeader;
  const std::vector< std::string > bksRow = lookupRow(dataDir + "/bks.tsv", instanceName, &bksHeader);
  if (!bksRow.empty()) {
    int bksCol = -1;
    for (size_t i = 0; i < bksHeader.size(); ++i) {
      if (bksHeader[i] == "bks" || bksHeader[i] == "objective") {
        bksCol = (int)i;
        break;
      }
    }
    if (bksCol < 0 && bksRow.size() > 1)
      bksCol = 1; // fall back to the simple 2-column schema's assumption
    if (bksCol >= 0 && bksCol < (int)bksRow.size()) {
      const double bks = atof(bksRow[bksCol].c_str());
      printf("[replay] %s: bks=%.10g\n", name.c_str(), bks);
    }
  }

  return 0;
}
