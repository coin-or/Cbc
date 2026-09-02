/**
 * Correctness/regression smoke test for the real CglImpliedClique class
 * (Cgl/src/CglImpliedClique), replaying it against the same fixture files
 * used by impclique-bench.cpp / impclique-vs-clique-bench.cpp. Loads a
 * fixture, calls CglImpliedClique::generateCuts once, applies any cuts
 * found, resolves, and reports the LP objective improvement -- so a
 * fixture already characterised by the ad hoc separate()/toRowCut() logic
 * in impclique-bench.cpp can be cross-checked against the production class
 * to confirm both produce the same cuts / bound movement.
 *
 * Usage:
 *   impclique-real-bench <stem> [--expect-improve=X] [--quiet]
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.cgraph, <stem>.bas.
 */

#include "CglImpliedClique.hpp"
#include "ClpSimplex.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace {

bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
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
      fprintf(stderr, "WARNING: %s: meta says padded but rows mismatch; leaving matrix alone\n",
        baseName(stem).c_str());
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
      fprintf(stderr, "WARNING: %s: no .ctype sidecar\n", baseName(stem).c_str());
    return -1;
  }
  int sidecarCols = -1;
  if (fscanf(fp, "cols %d\n", &sidecarCols) != 1 || sidecarCols != si.getNumCols()) {
    fprintf(stderr, "ERROR: %s: .ctype column-count mismatch\n", baseName(stem).c_str());
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

bool loadFixture(OsiClpSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";
  const std::string cgr = stem + ".cgraph";

  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no problem file for stem %s\n", stem.c_str());
    return false;
  }

  ClpSimplex *lp = si.getModelPtr();
  lp->setLogLevel(0);
  si.messageHandler()->setLogLevel(0);

  fflush(stdout);
  int savedStdout = dup(1);
  int devNull = open("/dev/null", O_WRONLY);
  if (devNull >= 0)
    dup2(devNull, 1);
  int mpsRc = si.readMps(mps.c_str());
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

  dropPadRow(si, stem, quiet);
  restoreColTypes(si, stem, quiet);

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
  if (cg->size() != (size_t)si.getNumCols() * 2) {
    fprintf(stderr, "ERROR: graph/model mismatch for %s\n", stem.c_str());
    delete cg;
    return false;
  }
  si.setCGraph(cg);

  return true;
}

std::string fixtureStem(const char *arg)
{
  std::string s(arg);
  static const char *suffixes[] = { ".mps.gz", ".mps", ".cgraph", ".bas", ".sol", ".bas.status" };
  for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
    const std::string suf(suffixes[i]);
    if (s.size() > suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0)
      return s.substr(0, s.size() - suf.size());
  }
  return s;
}

} // namespace

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <fixture-stem> [--expect-improve=X] [--quiet]\n", argv[0]);
    return 2;
  }

  bool quiet = false;
  double expectImprove = -1.0;
  double minViol = 0.02;
  std::string stemArg;
  for (int i = 1; i < argc; i++) {
    const std::string arg(argv[i]);
    if (arg == "--quiet") {
      quiet = true;
    } else if (arg.rfind("--expect-improve=", 0) == 0) {
      expectImprove = atof(arg.substr(strlen("--expect-improve=")).c_str());
    } else if (arg.rfind("--min-viol=", 0) == 0) {
      minViol = atof(arg.substr(strlen("--min-viol=")).c_str());
    } else if (arg.rfind("--", 0) != 0) {
      stemArg = arg;
    }
  }
  if (stemArg.empty()) {
    fprintf(stderr, "ERROR: missing fixture stem\n");
    return 2;
  }
  const std::string stem = fixtureStem(stemArg.c_str());

  OsiClpSolverInterface si;
  if (!loadFixture(si, stem, quiet))
    return 1;

  const double objStart = si.getObjValue();

  CglImpliedClique gen;
  gen.setMinViol(minViol);
  OsiCuts cs;
  gen.generateCuts(si, cs);

  const int nCuts = cs.sizeRowCuts();
  if (nCuts > 0) {
    si.applyCuts(cs);
    si.resolve();
    if (!si.isProvenOptimal()) {
      fprintf(stderr, "ERROR: %s: LP infeasible/not optimal after applying cuts\n", stem.c_str());
      return 1;
    }
  }
  const double objEnd = si.getObjValue();
  const double objImprove = objEnd - objStart;

  printf("%s,cuts=%d,objStart=%.6f,objEnd=%.6f,objImprove=%.6f\n",
    baseName(stem).c_str(), nCuts, objStart, objEnd, objImprove);

  if (expectImprove >= 0.0) {
    const double tol = 1e-3 * (fabs(expectImprove) > 1.0 ? fabs(expectImprove) : 1.0);
    if (fabs(objImprove - expectImprove) > tol) {
      fprintf(stderr, "FAIL: %s: expected objImprove ~%.6f, got %.6f\n",
        stem.c_str(), expectImprove, objImprove);
      return 1;
    }
  }

  return 0;
}
