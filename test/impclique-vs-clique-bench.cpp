/**
 * Marginal-bound benchmark: how much does the implication-clique separator
 * add *on top of* standard clique separation (CglBKClique), not in isolation?
 *
 * @file impclique-vs-clique-bench.cpp
 * @brief chain CglBKClique rounds, then this separator, on the same fixture
 *
 * Every cut this separator produces is literally a clique inequality in the
 * conflict graph (see impclique-explore.cpp's header comment), so the
 * natural question is whether a generic clique separator already finds them
 * once it is given enough rounds. CglBKClique roots its search at fractional
 * columns; this separator specifically roots at complement(y) for every
 * disaggregation-row hub y, fractional or not (a 0/1-valued y is exactly the
 * case CglBKClique's fractional-candidate pool skips). The two are not
 * mutually exclusive, so this tool measures them in sequence on one fixture:
 *
 *   objStart --[CglBKClique, N rounds]--> objAfterCliques
 *            --[this separator, mixed-literals]--> objFinal
 *
 * reporting both the standard-clique-only improvement and this separator's
 * marginal improvement on top of it -- the number that matters for deciding
 * whether to add this as a *new* generator versus expecting CglBKClique to
 * already cover it if just tuned differently.
 *
 * Fixture loading, CglBKClique round loop and defaults are copied verbatim
 * from bkclique-bench.cpp; separate()/toRowCut are copied verbatim from
 * impclique-bench.cpp (mixed-literals variant). Kept as one flat file (like
 * the other *-bench tools) rather than factored into a shared header, since
 * these benches are meant to be readable standalone and are not part of the
 * production build.
 *
 * Usage:
 *   impclique-vs-clique-bench <stem> [--rounds=N] [--csv-header] [--quiet]
 *   impclique-vs-clique-bench --header
 */

#include "CglBKClique.hpp"
#include "CglTreeInfo.hpp"
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
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

double wallClock() { return CoinGetTimeOfDay(); }

bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

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
  si.getColType(true);
  return restored;
}

struct Fixture {
  OsiClpSolverInterface si;
  bool ok = false;
};

bool loadFixture(Fixture &f, const std::string &stem, bool quiet)
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
  restoreColTypes(f.si, stem, quiet);

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

  if (haveBasis) {
    f.si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
    f.si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
    f.si.resolve();
  } else {
    f.si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);
    f.si.initialSolve();
  }

  if (!f.si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: LP not optimal after warm start (%s)\n", stem.c_str());
    return false;
  }

  if (!fileExists(cgr)) {
    fprintf(stderr, "ERROR: no conflict graph %s\n", cgr.c_str());
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

  f.ok = true;
  return true;
}

const double EPS = 1e-6;

/// A clique-cut member is a literal, not just a column -- see
/// impclique-bench.cpp's Lit for the full rationale (mixed-literals mode).
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
};

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

/// Implication-clique separator, mixed-literals variant (see
/// impclique-bench.cpp for the identical logic and rationale).
std::vector< FoundCut > separate(const OsiSolverInterface &si)
{
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

    auto conf = cg->conflictingNodes(yComp, temp.data(), iv.data());

    std::vector< std::pair< double, int > > candidates;
    for (size_t k = 0; k < conf.first; k++) {
      const size_t node = conf.second[k];
      int col;
      bool neg;
      if (node < static_cast< size_t >(numCols)) {
        col = static_cast< int >(node);
        neg = false;
      } else if (node < static_cast< size_t >(2 * numCols)) {
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
    if (candidates.empty())
      continue;

    std::sort(candidates.begin(), candidates.end(),
      [](const std::pair< double, int > &a, const std::pair< double, int > &b) {
        return a.first > b.first;
      });

    std::vector< int > cliqueNodes = greedyClique(cg, candidates);

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

void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "       %s --header\n"
    "\n"
    "Runs two independent round-by-round trajectories from the same starting\n"
    "LP, mirroring how CBC's cut loop actually works (every active generator\n"
    "fires against the current LP each pass, then one resolve, repeat):\n"
    "  A (clique-only):     each round, CglBKClique only\n"
    "  B (clique+implied):  each round, CglBKClique AND this separator\n"
    "                       together, their cuts applied in the same resolve\n"
    "A single chain-after-convergence run (cliques to exhaustion, then this\n"
    "separator once) would overstate this separator's marginal value, since\n"
    "a cut it finds in round 1 can change what CglBKClique finds in round 2\n"
    "and vice versa -- this tool measures the interleaved, realistic case.\n"
    "\n"
    "Options:\n"
    "  --rounds=N    max rounds for each trajectory (default 4)\n"
    "  --csv-header  print the CSV header before the data line\n"
    "  --quiet       suppress warnings\n",
    prog, prog);
}

const char *CSV_HEADER
  = "name,maxRounds,roundsA,roundsB,cliqueCutsA,cliqueCutsB,impCutsB,"
    "objStart,objA,objB,cliqueImprove,combinedImprove,marginalImprove,"
    "marginalFracOfCombined,perRoundObjA,perRoundObjB";

CglBKClique makeBkClique()
{
  CglBKClique bk;
  bk.setMaxCallsBK(1000);
  bk.setExtendingMethod(4);
  bk.setPivotingStrategy(CoinBronKerbosch::PivotingStrategy::Weight);
  bk.setMaxInducedSize(10000);
  bk.setMinViol(0.02);
  return bk;
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

  bool csvHeader = false, quiet = false;
  int maxRounds = 4;
  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (!strncmp(a, "--rounds=", 9))
      maxRounds = atoi(a + 9);
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

  // Two independent fixture loads -- deliberately not a clone of one solver,
  // to avoid any question of shared CGraph ownership between the two
  // trajectories once they diverge after round 1.
  Fixture fA, fB;
  if (!loadFixture(fA, stem, quiet) || !loadFixture(fB, stem, quiet))
    return 1;

  const double objStart = fA.si.getObjValue();
  const double sense = fA.si.getObjSenseInCbc();

  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  info.formulation_rows = fA.si.getNumRows();
  info.inTree = false;
  info.options = 0;

  bool aAlive = true, bAlive = true;
  double objA = objStart, objB = objStart;
  int roundsA = 0, roundsB = 0;
  int cliqueCutsA = 0, cliqueCutsB = 0, impCutsB = 0;
  std::string perRoundA, perRoundB;

  for (int round = 0; round < maxRounds && (aAlive || bAlive); ++round) {
    if (aAlive) {
      CglBKClique bk = makeBkClique();
      OsiCuts cs;
      bk.generateCuts(fA.si, cs, info);
      const int nCuts = cs.sizeRowCuts();
      if (nCuts == 0) {
        aAlive = false;
      } else {
        cliqueCutsA += nCuts;
        fA.si.applyCuts(cs);
        fA.si.resolve();
        if (fA.si.isProvenOptimal()) {
          const double objRoundEnd = fA.si.getObjValue();
          char buf[64];
          snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "",
            sense * (objRoundEnd - objA));
          perRoundA += buf;
          objA = objRoundEnd;
          roundsA++;
        } else {
          aAlive = false;
        }
      }
    }

    if (bAlive) {
      // Both generators fire against the *same* current LP, exactly as
      // CBC's cut loop calls every active generator once per pass -- neither
      // sees the other's cuts until after the shared resolve.
      CglBKClique bk = makeBkClique();
      OsiCuts cs;
      bk.generateCuts(fB.si, cs, info);
      const std::vector< FoundCut > impCuts = separate(fB.si);

      const int nClique = cs.sizeRowCuts();
      const int nImp = static_cast< int >(impCuts.size());
      if (nClique == 0 && nImp == 0) {
        bAlive = false;
      } else {
        cliqueCutsB += nClique;
        impCutsB += nImp;
        for (const auto &c : impCuts)
          cs.insert(toRowCut(c));
        fB.si.applyCuts(cs);
        fB.si.resolve();
        if (fB.si.isProvenOptimal()) {
          const double objRoundEnd = fB.si.getObjValue();
          char buf[64];
          snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "",
            sense * (objRoundEnd - objB));
          perRoundB += buf;
          objB = objRoundEnd;
          roundsB++;
        } else {
          bAlive = false;
        }
      }
    }
  }

  const double cliqueImprove = sense * (objA - objStart);
  const double combinedImprove = sense * (objB - objStart);
  const double marginalImprove = sense * (objB - objA);
  const double marginalFrac
    = fabs(combinedImprove) > 1e-9 ? marginalImprove / combinedImprove : 0.0;

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%d,%d,%d,%d,%.15g,%.15g,%.15g,%.6g,%.6g,%.6g,%.6g,%s,%s\n",
    baseName(stem).c_str(), maxRounds, roundsA, roundsB, cliqueCutsA,
    cliqueCutsB, impCutsB, objStart, objA, objB, cliqueImprove,
    combinedImprove, marginalImprove, marginalFrac,
    perRoundA.empty() ? "0" : perRoundA.c_str(),
    perRoundB.empty() ? "0" : perRoundB.c_str());

  return 0;
}
