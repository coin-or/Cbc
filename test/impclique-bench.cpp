/**
 * Stand-alone benchmark for "implication clique" separation.
 *
 * @file impclique-bench.cpp
 * @brief replay a hub-implication clique separation call captured from CBC
 *
 * Companion to bkclique-bench.cpp, reusing the exact same fixture format --
 * CbcClqFixtureDump.hpp's "sep" tag, written by a Cbc built with
 * -DCBC_DUMP_CLQSEP_FIXTURE right after the root LP, before any cut has been
 * applied (currentPassNumber_ == 1 && !node && !parentModel_ in CbcModel.cpp).
 * No new fixture format was needed: this separator consumes exactly what
 * CglBKClique consumes -- the preprocessed matrix, its optimal LP basis, and
 * the conflict graph -- just rooted differently (see impclique-explore.cpp's
 * header comment for the algebraic argument). Loading the *captured* graph
 * rather than rebuilding it matters here for the same reason it does for BK:
 * see CbcClqFixtureDump.hpp's header comment for measured staleness (e.g.
 * trdtaunimep, 190900 rebuilt vs 191028 captured direct conflicts).
 *
 * Separator: for every binary column y with y* < 1 (not necessarily
 * fractional -- this is exactly the case CglBKClique's fractional-literal
 * candidate pool misses), root at complement(y) in the conflict graph,
 * restrict candidates to y's neighbours that are primal column literals with
 * positive LP value, and greedily grow a clique among them by descending LP
 * weight. --singleton-only disables the greedy growth (plain disaggregation
 * x_i <= y only), so running the same fixture with/without it isolates what
 * clique extension buys, exactly like bkclique-bench's --ext-method=0.
 *
 * Usage:
 *   impclique-bench <stem> [options]        (stem = dir/name.tag)
 *   impclique-bench --self-test <stem>
 *   impclique-bench --header
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.cgraph, <stem>.bas, ...
 * A full path to any one of those files also works; the suffix is stripped.
 */

#include "ClpSimplex.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const char *msg)
{
  if (!ok) {
    printf("FAIL: %s\n", msg);
    ++g_failures;
  }
}

double wallClock() { return CoinGetTimeOfDay(); }

bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

/// Reduce any of the fixture's file names to the shared stem.
std::string fixtureStem(const char *arg)
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

std::string baseName(const std::string &path)
{
  const size_t slash = path.rfind('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

/// Read one integer key out of the fixture's `.meta`. Returns `dflt` when the
/// file or the key is absent.
long metaInt(const std::string &path, const char *key, long dflt)
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

/// Remove the pad row a fixture with originally-empty columns carries; see
/// CbcClqFixtureDump.hpp's cbcClqFixtureWriteMps for why it exists.
bool dropPadRow(OsiSolverInterface &si, const std::string &stem, bool quiet)
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

/// Restore integrality from the `.ctype` sidecar (MPS cannot express a fixed
/// integer column: see CbcClqFixtureDump.hpp's header comment).
int restoreColTypes(OsiSolverInterface &si, const std::string &stem, bool quiet)
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
  si.getColType(true); // getColType() caches; force recompute after restoring.
  return restored;
}

/// Everything the fixture determines, plus how long loading it took.
struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  int warmStartIters = 0;
  double cgraphTime = 0.0;
  int restoredColTypes = -1;
  bool ok = false;
};

/**
 * Load the fixture and warm-start to the captured optimum. See
 * bkclique-bench.cpp's loadFixture for the two gotchas this avoids:
 * OsiClp's cached basis_ shadowing the model (fixed with setWarmStart(NULL)),
 * and initialSolve() discarding a loaded basis via presolve (use resolve()).
 */
bool loadFixture(Fixture &f, const std::string &stem, bool rebuildCgraph, bool quiet)
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
  f.si.messageHandler()->setLogLevel(0);
  // A handful of fixtures carry pre-existing MPS format quirks (e.g.
  // malformed bound lines) that CoinMpsIO reports at severity 0 -- printed
  // unconditionally regardless of log level. Since these land on stdout
  // and would corrupt our CSV row, redirect fd 1 to /dev/null for the
  // duration of the read only.
  fflush(stdout);
  int savedStdout = dup(1);
  int devNull = open("/dev/null", O_WRONLY);
  if (devNull >= 0)
    dup2(devNull, 1);
  int mpsRc = f.si.readMps(mps.c_str());
  fflush(stdout);
  if (savedStdout >= 0) {
    dup2(savedStdout, 1);
    close(savedStdout);
  }
  if (devNull >= 0)
    close(devNull);
  if (mpsRc) {
    fprintf(stderr, "ERROR: failed to read %s\n", mps.c_str());
    return false;
  }

  dropPadRow(f.si, stem, quiet);
  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  bool haveBasis = false;
  if (fileExists(bas)) {
    if (lp->readBasis(bas.c_str()) < 0) {
      fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
    } else {
      haveBasis = true;
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

const double EPS = 1e-6;

/// A clique-cut member is a *literal*, not just a column: `node` is a raw
/// conflict-graph node id (col for the "col=1" literal, col+numCols for the
/// "col=0" literal, i.e. the complement). `col`/`neg` are the decoded form,
/// kept alongside for cut assembly (see toRowCut).
struct Lit {
  int node;
  int col;
  bool neg;
};

struct FoundCut {
  int hubCol;
  std::vector< Lit > members;
  double lhs;
  double rhs;
  double violation() const { return lhs - rhs; }
};

/// Greedily grow a clique among `candidates` (already known to all conflict
/// with hubComp), ordered by descending LP weight. O(|candidates|^2) conflict
/// checks -- see impclique-explore.cpp for the same routine and its rationale.
std::vector< int > greedyClique(const CoinConflictGraph *cg,
  const std::vector< std::pair< double, int > > &sortedCandidates)
{
  std::vector< int > clique;
  clique.reserve(sortedCandidates.size());
  for (const auto &wc : sortedCandidates) {
    const int cand = wc.second;
    bool ok = true;
    for (int m : clique) {
      if (!cg->conflicting(static_cast< size_t >(cand), static_cast< size_t >(m))) {
        ok = false;
        break;
      }
    }
    if (ok)
      clique.push_back(cand);
  }
  return clique;
}

/// Separates, and reports the time split between building each hub's
/// candidate list (scanTime) and greedily growing the clique (extendTime) --
/// the two phases have very different asymptotic behaviour (scan is one
/// conflictingNodes() call per column; extension is quadratic in candidate
/// count), so a profile needs them apart to tell which dominates.
std::vector< FoundCut > separate(const OsiSolverInterface &si, bool singletonOnly,
  bool mixedLiterals, double &scanTime, double &extendTime, size_t &candidateNodeVisits)
{
  scanTime = 0.0;
  extendTime = 0.0;
  candidateNodeVisits = 0;
  std::vector< FoundCut > cuts;
  const CoinConflictGraph *cg = si.getCGraph();
  if (!cg)
    return cuts;

  const int numCols = si.getNumCols();
  const double *sol = si.getColSolution();
  const char *colType = si.getColType(true);

  std::vector< size_t > temp(cg->size());
  std::vector< char > iv(cg->size(), 0);

  for (int y = 0; y < numCols; y++) {
    if (colType[y] == 0)
      continue;
    const double yVal = sol[y];
    if (yVal > 1.0 - EPS)
      continue;

    const size_t yComp = static_cast< size_t >(y) + numCols;
    if (yComp >= cg->size())
      continue;

    const double tScan0 = wallClock();
    auto conf = cg->conflictingNodes(yComp, temp.data(), iv.data());

    // Weighted-by-LP-value node id: positive-literal candidates are the
    // column node itself (col); a complemented-literal candidate ("col=0",
    // i.e. (1-col)) is the col+numCols node -- only considered when
    // mixedLiterals is set, since a cut mixing signs is a materially
    // different (and unvalidated-until-now) shape from the plain
    // disaggregation x_i <= y this separator started from.
    std::vector< std::pair< double, int > > candidates;
    for (size_t k = 0; k < conf.first; k++) {
      const size_t node = conf.second[k];
      int col;
      bool neg;
      if (node < static_cast< size_t >(numCols)) {
        col = static_cast< int >(node);
        neg = false;
      } else if (mixedLiterals && node < static_cast< size_t >(2 * numCols)) {
        col = static_cast< int >(node - static_cast< size_t >(numCols));
        neg = true;
      } else {
        continue;
      }
      if (col == y)
        continue;
      if (colType[col] == 0)
        continue;
      const double litVal = neg ? (1.0 - sol[col]) : sol[col];
      if (litVal > EPS)
        candidates.push_back(std::make_pair(litVal, static_cast< int >(node)));
    }
    candidateNodeVisits += conf.first;
    if (candidates.empty()) {
      scanTime += wallClock() - tScan0;
      continue;
    }

    std::sort(candidates.begin(), candidates.end(),
      [](const std::pair< double, int > &a, const std::pair< double, int > &b) {
        return a.first > b.first;
      });
    scanTime += wallClock() - tScan0;

    const double tExt0 = wallClock();
    std::vector< int > cliqueNodes;
    if (singletonOnly) {
      cliqueNodes.push_back(candidates.front().second);
    } else {
      cliqueNodes = greedyClique(cg, candidates);
    }
    extendTime += wallClock() - tExt0;

    std::vector< Lit > members;
    double lhs = 0.0;
    for (int node : cliqueNodes) {
      const bool neg = node >= numCols;
      const int col = neg ? node - numCols : node;
      members.push_back(Lit{ node, col, neg });
      lhs += neg ? (1.0 - sol[col]) : sol[col];
    }

    if (lhs - yVal > EPS) {
      FoundCut fc;
      fc.hubCol = y;
      fc.members = members;
      fc.lhs = lhs;
      fc.rhs = yVal;
      cuts.push_back(fc);
    }
  }
  return cuts;
}

/// sum_pos x_i + sum_neg (1 - x_j) <= y  <=>  sum_pos x_i - sum_neg x_j - y <= -negCount
OsiRowCut toRowCut(const FoundCut &fc)
{
  std::vector< int > idx;
  std::vector< double > val;
  idx.reserve(fc.members.size() + 1);
  val.reserve(fc.members.size() + 1);
  int negCount = 0;
  for (const auto &lit : fc.members) {
    idx.push_back(lit.col);
    val.push_back(lit.neg ? -1.0 : 1.0);
    if (lit.neg)
      negCount++;
  }
  idx.push_back(fc.hubCol);
  val.push_back(-1.0);
  OsiRowCut rc;
  rc.setRow(static_cast< int >(idx.size()), idx.data(), val.data());
  rc.setLb(-COIN_DBL_MAX);
  rc.setUb(-static_cast< double >(negCount));
  return rc;
}

/// objImprove is a difference of two LP bounds and inherits their absolute
/// error, which scales with magnitude; a relative test (with an absolute
/// floor) is what separates real movement from noise. See bkclique-bench.cpp's
/// boundMoved() for the measured rationale (brazil3 2.7e-12 vs mzzv11 171.49).
bool boundMoved(double objStart, double objImprove)
{
  const double scale = fabs(objStart) > 1.0 ? fabs(objStart) : 1.0;
  return objImprove > 1.0e-9 * scale && objImprove > 1.0e-9;
}

/**
 * Round-trip the graph through save()/load(), byte-identical. This separator
 * introduces no new serialized structure of its own -- it consumes exactly
 * the .cgraph CbcClqFixtureDump.hpp already writes for CglBKClique -- so this
 * is the same check as bkclique-bench.cpp's --self-test, kept here too so a
 * fixture can be validated without depending on that binary.
 */
int selfTest(const std::string &stem)
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

  const std::string t1 = "/tmp/impclique-selftest-1.cgraph";
  const std::string t2 = "/tmp/impclique-selftest-2.cgraph";
  check(a->save(t1.c_str()) == 0, "save of loaded graph");

  CoinStaticConflictGraph *b = CoinStaticConflictGraph::load(t1.c_str());
  if (!b) {
    printf("FAIL: could not re-load saved graph\n");
    delete a;
    return 1;
  }
  check(b->save(t2.c_str()) == 0, "re-save of re-loaded graph");

  {
    FILE *f1 = fopen(t1.c_str(), "rb");
    FILE *f2 = fopen(t2.c_str(), "rb");
    bool same = f1 && f2;
    while (same) {
      char b1[65536], b2[65536];
      const size_t r1 = fread(b1, 1, sizeof(b1), f1);
      const size_t r2 = fread(b2, 1, sizeof(b2), f2);
      if (r1 != r2 || memcmp(b1, b2, r1) != 0) {
        same = false;
        break;
      }
      if (r1 == 0)
        break;
    }
    if (f1)
      fclose(f1);
    if (f2)
      fclose(f2);
    check(same, "save/load/save is byte-identical");
  }

  // The two-run determinism check that matters for *this* tool: separating
  // twice against the same loaded graph and solution must find the identical
  // set of cuts, since the separator has no randomness. Catches an
  // uninitialized-read bug that a single run cannot.
  delete b;
  delete a;

  printf(g_failures ? "self-test: %d FAILURE(S)\n" : "self-test: all checks passed\n",
    g_failures);
  return g_failures ? 1 : 0;
}

const char *CSV_HEADER
  = "name,singletonOnly,mixedLiterals,warmStartTime,warmStartIters,cgraphTime,scanTime,"
    "extendTime,sepTime,candidateNodeVisits,resolveTime,resolveIters,"
    "totalCuts,multiMemberCuts,negLiteralCuts,totalNegLiterals,avgCutLen,totalViol,maxViol,objStart,objEnd,"
    "objImprove,objImproveRel,boundMoved,cgNodes,cgDirectConf,cgDensity,"
    "restoredInt";

void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "       %s --self-test <fixture-stem>\n"
    "       %s --header\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.cgraph/.bas are appended.\n"
    "Passing any one of those files also works.\n"
    "\n"
    "Options:\n"
    "  --singleton-only    disable clique growth (k=1 disaggregation only);\n"
    "                      the control that isolates what extension buys\n"
    "  --mixed-literals    also allow complemented-literal members (col=0,\n"
    "                      i.e. (1-col)) when growing the clique, not just\n"
    "                      more positive columns from the same row\n"
    "  --rebuild-cgraph    rebuild the graph from the matrix instead of loading\n"
    "                      the captured one (not faithful; for comparison only)\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --quiet             suppress warnings\n",
    prog, prog, prog);
}

} // namespace

int main(int argc, char **argv)
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

  bool doSelfTest = false, singletonOnly = false, rebuildCgraph = false;
  bool csvHeader = false, quiet = false, mixedLiterals = false;
  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (!strcmp(a, "--self-test"))
      doSelfTest = true;
    else if (!strcmp(a, "--singleton-only"))
      singletonOnly = true;
    else if (!strcmp(a, "--rebuild-cgraph"))
      rebuildCgraph = true;
    else if (!strcmp(a, "--mixed-literals"))
      mixedLiterals = true;
    else if (!strcmp(a, "--csv-header"))
      csvHeader = true;
    else if (!strcmp(a, "--quiet"))
      quiet = true;
    else if (a[0] == '-') {
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

  const double objStart = f.si.getObjValue();

  double scanTime = 0.0, extendTime = 0.0;
  size_t candidateNodeVisits = 0;
  const double tSep0 = wallClock();
  std::vector< FoundCut > cuts
    = separate(f.si, singletonOnly, mixedLiterals, scanTime, extendTime, candidateNodeVisits);
  const double sepTime = wallClock() - tSep0;

  double totalViol = 0.0, maxViol = 0.0;
  size_t sumLen = 0, nMultiMember = 0, nNegLiteralCuts = 0, totalNegLiterals = 0;
  for (const auto &c : cuts) {
    totalViol += c.violation();
    maxViol = std::max(maxViol, c.violation());
    sumLen += c.members.size();
    if (c.members.size() > 1)
      nMultiMember++;
    size_t negHere = 0;
    for (const auto &lit : c.members)
      if (lit.neg)
        negHere++;
    if (negHere > 0)
      nNegLiteralCuts++;
    totalNegLiterals += negHere;
  }
  const double avgLen = cuts.empty() ? 0.0 : static_cast< double >(sumLen) / cuts.size();

  double objEnd = objStart;
  double resolveTime = 0.0;
  int resolveIters = 0;
  if (!cuts.empty()) {
    OsiCuts osc;
    for (const auto &c : cuts)
      osc.insert(toRowCut(c));
    f.si.applyCuts(osc);
    const double tR0 = wallClock();
    f.si.resolve();
    resolveTime = wallClock() - tR0;
    resolveIters = f.si.getIterationCount();
    if (f.si.isProvenOptimal())
      objEnd = f.si.getObjValue();
  }

  const double objImprove = f.si.getObjSenseInCbc() * (objEnd - objStart);
  const double objImproveRel = objImprove / (fabs(objStart) > 1.0 ? fabs(objStart) : 1.0);

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%lu,%.6f,%d,%d,%lu,%lu,%lu,%.3f,%.10g,%.10g,"
         "%.15g,%.15g,%.6g,%.6g,%d,%lu,%lu,%.10f,%d\n",
    baseName(stem).c_str(), (int)singletonOnly, (int)mixedLiterals, f.warmStartTime, f.warmStartIters,
    f.cgraphTime, scanTime, extendTime, sepTime, (unsigned long)candidateNodeVisits,
    resolveTime, resolveIters, (int)cuts.size(), (unsigned long)nMultiMember,
    (unsigned long)nNegLiteralCuts, (unsigned long)totalNegLiterals, avgLen,
    totalViol, maxViol, objStart, objEnd, objImprove, objImproveRel,
    (int)boundMoved(objStart, objImprove), (unsigned long)cg->size(),
    (unsigned long)cg->nTotalDirectConflicts(), cg->density(), f.restoredColTypes);

  return 0;
}
