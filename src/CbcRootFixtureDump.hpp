/**
 * Fixture dumping for root-node processing experiments (cut generation +
 * heuristics), generator- and heuristic-agnostic.
 *
 * @file CbcRootFixtureDump.hpp
 * @brief dump the preprocessed problem + optimal root LP basis, for offline
 *        replay of the whole root-node cut generation / heuristics loop
 *
 * Root-node processing -- pre-processing, then the cut generation loop, then
 * root heuristics -- only becomes interesting to experiment on *after*
 * pre-processing and the root LP solve have both happened: that is the state
 * every cut generator and every root heuristic actually starts from. Reaching
 * that state costs a full pre-processing pass + root LP solve per experiment,
 * which on some mip-sanity-data instances is itself the expensive part. This
 * writes out everything needed to reconstruct that state directly, so a replay
 * driver can rebuild a CbcModel around it and re-enter the same root loop
 * (solveWithCuts + doHeuristicsAtRoot) without paying for pre-processing or the
 * root LP again:
 *
 *   <name>.<tag>.mps.gz    the preprocessed problem, formatType 2 (IEEE hex)
 *   <name>.<tag>.bas       the optimal root LP basis (warm start, 0 iterations)
 *   <name>.<tag>.sol       the optimal root LP solution, for cross-checking
 *   <name>.<tag>.ctype     integer/continuous per column (MPS cannot say: see
 *                          CbcClqFixtureDump.hpp for why this matters)
 *   <name>.<tag>.meta      rows/cols/density/objValue/lpOptimal/... so a driver
 *                          can filter fixtures without loading them
 *
 * Unlike the per-generator fixtures (CbcClqFixtureDump.hpp and friends), this
 * one carries no generator-specific payload at all -- a replay driver is
 * expected to attach whatever cut generators / heuristics it wants to
 * experiment with itself, exactly as the normal `cbc` driver would for a fresh
 * run, rather than replaying a single call captured mid-loop. That is the
 * point: one fixture per instance is enough to explore any cut generator /
 * heuristic strategy, rather than one fixture set per generator.
 *
 * Captured at the same call site as the per-generator fixtures, for the same
 * reason: right after the root LP has been solved and before any cut has been
 * installed, which is the one point where the preprocessed matrix and the LP
 * solution/basis are simultaneously consistent and nothing generator-specific
 * has run yet.
 *
 * Entirely behind CBC_DUMP_ROOT_FIXTURE and off by default: build with
 * -DCBC_DUMP_ROOT_FIXTURE to generate fixtures. Header-only static functions,
 * so no Makefile.am/Makefile.in changes are needed to carry it.
 *
 * Environment:
 *   CBC_ROOT_FIXTURE_DIR   output directory
 *                          (default ~/instances/mip-sanity-data/rootFixtures)
 *   CBC_ROOT_FIXTURE_NAME  instance base name, when the MPS problem name is
 *                          unhelpful (drivers normally set this)
 **/

#ifndef CbcRootFixtureDump_H
#define CbcRootFixtureDump_H

#ifdef CBC_DUMP_ROOT_FIXTURE

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

/// mkdir -p, so a driver need not pre-create the tree.
static void cbcRootFixtureMkdirP(const std::string &path)
{
  for (size_t i = 1; i <= path.size(); ++i) {
    if (i < path.size() && path[i] != '/')
      continue;
    const std::string part = path.substr(0, i);
#ifdef _WIN32
    _mkdir(part.c_str());
#else
    mkdir(part.c_str(), 0775);
#endif
  }
}

/// Output directory, from the environment or the default fixture location.
static std::string cbcRootFixtureDir()
{
  const char *env = getenv("CBC_ROOT_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/mip-sanity-data/rootFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback,
/// but drivers should set CBC_ROOT_FIXTURE_NAME: several instances share a
/// problem name, and some carry none at all.
static std::string cbcRootFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_ROOT_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

/// Collect the model's row and column names into the `const char **` pair
/// writeMpsNative wants, so the MPS carries the same names the `.bas` will.
/// getRowName/getColName are non-const, hence the non-const solver argument.
static void cbcRootFixtureNames(OsiSolverInterface *si,
  std::vector< std::string > &rowStore, std::vector< std::string > &colStore,
  std::vector< const char * > &rowPtrs, std::vector< const char * > &colPtrs)
{
  const int nRows = si->getNumRows(), nCols = si->getNumCols();
  rowStore.resize(nRows);
  colStore.resize(nCols);
  rowPtrs.resize(nRows);
  colPtrs.resize(nCols);
  for (int i = 0; i < nRows; ++i) {
    rowStore[i] = si->getRowName(i);
    rowPtrs[i] = rowStore[i].c_str();
  }
  for (int j = 0; j < nCols; ++j) {
    colStore[j] = si->getColName(j);
    colPtrs[j] = colStore[j].c_str();
  }
}

/**
 * Write the problem, preserving column *indices* (a replay driver reapplies
 * `.ctype`/`.bas` by index, so a shifted column would silently mismatch them).
 *
 * Same empty-column padding trick as CbcClqFixtureDump.hpp: CoinMpsIO::writeMps
 * only emits a column that has matrix elements or a nonzero objective, so an
 * empty column would otherwise vanish and shift every later index down. See
 * that file's header comment for the measured impact and the round-trip check;
 * the fix here is identical, just without the conflict-graph angle.
 */
static int cbcRootFixtureWriteMps(OsiSolverInterface *si, const std::string &path,
  int &paddedColumns)
{
  const int nCols = si->getNumCols();
  std::vector< int > empties;

  const CoinPackedMatrix *byCol = si->getMatrixByCol();
  if (byCol && si->getNumRows() > 0) {
    const int *lengths = byCol->getVectorLengths();
    const double *obj = si->getObjCoefficients();
    for (int j = 0; j < nCols; ++j) {
      if (lengths[j] == 0 && obj[j] == 0.0)
        empties.push_back(j);
    }
  }
  paddedColumns = (int)empties.size();

  std::vector< std::string > rowStore, colStore;
  std::vector< const char * > rowPtrs, colPtrs;

  if (empties.empty()) {
    // formatType 2 = IEEE hex, so every coefficient round-trips bit-for-bit;
    // CoinMpsIO adds free format itself when a name exceeds 8 characters.
    cbcRootFixtureNames(si, rowStore, colStore, rowPtrs, colPtrs);
    return si->writeMpsNative(path.c_str(),
      rowPtrs.empty() ? NULL : &rowPtrs[0],
      colPtrs.empty() ? NULL : &colPtrs[0], 2, 1);
  }

  OsiSolverInterface *clone = si->clone();
  if (!clone)
    return -1;

  const double *cl = clone->getColLower();
  const double *cu = clone->getColUpper();
  double padLb = 0.0, padUb = 0.0;
  bool finite = true;
  for (size_t k = 0; k < empties.size(); ++k) {
    const double lo = cl[empties[k]], up = cu[empties[k]];
    if (lo <= -1.0e30 || up >= 1.0e30) {
      finite = false;
      break;
    }
    padLb += lo < 0.0 ? lo : 0.0;
    padUb += up > 0.0 ? up : 0.0;
  }
  if (!finite) {
    padLb = -1.0e29;
    padUb = 1.0e29;
  }

  const std::vector< double > coefs(empties.size(), 1.0);
  clone->addRow((int)empties.size(), &empties[0], &coefs[0], padLb, padUb);

  // Named after the pad row is added, so the extra row gets a name too; a
  // consumer deletes it anyway (see `.meta`'s `paddedColumns`).
  cbcRootFixtureNames(clone, rowStore, colStore, rowPtrs, colPtrs);
  const int rc = clone->writeMpsNative(path.c_str(),
    rowPtrs.empty() ? NULL : &rowPtrs[0],
    colPtrs.empty() ? NULL : &colPtrs[0], 2, 1);
  delete clone;
  return rc;
}

/**
 * Write the basis. OsiClp's writeBasisNative() emits FREEIEEE, which round-trips
 * doubles exactly. The base class implementation is a no-op that still returns
 * success, so success is confirmed by the file existing rather than by the
 * return code; failing that, fall back to a plain status dump good enough to
 * rebuild a CoinWarmStartBasis.
 */
static bool cbcRootFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
{
  si->writeBasisNative(path.c_str());
  {
    FILE *probe = fopen(path.c_str(), "r");
    if (probe) {
      fclose(probe);
      return true;
    }
  }

  const CoinWarmStartBasis *ws
    = dynamic_cast< const CoinWarmStartBasis * >(si->getWarmStart());
  if (!ws)
    return false;

  FILE *fp = fopen((path + ".status").c_str(), "w");
  if (!fp) {
    delete ws;
    return false;
  }
  fprintf(fp, "STATUS %d %d\n", ws->getNumStructural(), ws->getNumArtificial());
  for (int j = 0; j < ws->getNumStructural(); ++j)
    fprintf(fp, "C %d %d\n", j, (int)ws->getStructStatus(j));
  for (int i = 0; i < ws->getNumArtificial(); ++i)
    fprintf(fp, "R %d %d\n", i, (int)ws->getArtifStatus(i));
  fclose(fp);
  delete ws;
  return true;
}

/// Write the LP solution: an "=obj=" line, then one line per nonzero as
/// "<index> <name> <value>" -- same shape as the existing preProcessedInstances
/// fixture set.
static bool cbcRootFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  fprintf(fp, "=obj= %.15g\n", si->getObjValue());
  const double *x = si->getColSolution();
  if (x) {
    const int n = si->getNumCols();
    for (int j = 0; j < n; ++j) {
      if (x[j] == 0.0)
        continue;
      fprintf(fp, "%5d %-24s %.15g\n", j, si->getColName(j).c_str(), x[j]);
    }
  }
  fclose(fp);
  return true;
}

/**
 * Write the column types, which the MPS cannot carry for fixed integers.
 * See CbcClqFixtureDump.hpp's header comment for why this sidecar exists and
 * the measured impact of skipping it (39994/72141 integer columns lost on
 * physiciansched3-3, all of them fixed).
 */
static bool cbcRootFixtureWriteColTypes(const OsiSolverInterface *si,
  const std::string &path, int &integerColumns)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  const int n = si->getNumCols();
  const char *ct = si->getColType(true);
  integerColumns = 0;
  fprintf(fp, "cols %d\n", n);
  for (int j = 0; j < n; ++j) {
    if (si->isContinuous(j))
      continue;
    ++integerColumns;
    fprintf(fp, "%d %d\n", j, (int)ct[j]);
  }
  fclose(fp);
  return true;
}

/**
 * Write the provenance file. Carries only what a replay driver needs to sanity
 * check a fixture without loading it, plus the fields a rebuilt CbcModel cannot
 * otherwise infer (objective sense, whether the root LP was actually optimal).
 */
static bool cbcRootFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  int paddedColumns, int integerColumns, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  fprintf(fp, "tag %s\n", tag);
  fprintf(fp, "rows %d\n", si->getNumRows());
  fprintf(fp, "cols %d\n", si->getNumCols());
  fprintf(fp, "elements %d\n", si->getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si->isProvenOptimal());
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index; see
  // cbcRootFixtureWriteMps. `rows` above is the captured count, so a consumer
  // that loads rows+1 rows should delete the last row to recover the captured
  // model.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  // How many columns the .ctype sidecar lists. A loader that restores fewer
  // than this has a stale or truncated sidecar.
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());
  fclose(fp);
  return true;
}

/**
 * Dump one root fixture. Returns true if files were written.
 *
 * Unlike the per-generator fixtures there is no precondition to gate on -- a
 * root LP always exists at this call site (or the model is already integer
 * feasible with no LP to resolve, handled below) -- so nothing is skipped;
 * every instance in a driven set gets a fixture.
 *
 * \param si  solver whose preprocessed matrix and root LP solution/basis are
 *            captured
 * \param tag distinguishes capture sites; "root" for the normal one
 */
static bool cbcDumpRootFixture(OsiSolverInterface *si, const char *tag)
{
  if (!si)
    return false;

  const std::string name = cbcRootFixtureName(si);
  const std::string dir = cbcRootFixtureDir();
  cbcRootFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcRootFixtureWriteMps(si, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk = cbcRootFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const bool metaOk = cbcRootFixtureWriteMeta(si, tag, paddedColumns, integerColumns,
    stem + ".meta");
  // Only when there is an LP to capture. isProvenOptimal() is false, e.g., when
  // pre-processing itself fixed every variable and left nothing to solve.
  const bool haveLp = si->isProvenOptimal();
  const bool basOk = haveLp ? cbcRootFixtureWriteBasis(si, stem + ".bas") : true;
  const bool solOk = haveLp ? cbcRootFixtureWriteSol(si, stem + ".sol") : true;

  printf("[rootfixture] %s.%s: %s rows=%d cols=%d int=%d padded=%d lpOptimal=%d "
         "obj=%.15g (mps=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, paddedColumns,
    (int)haveLp, si->getObjValue(),
    mpsRc, (int)ctOk, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  return mpsRc == 0 && ctOk && metaOk && basOk && solOk;
}

#endif // CBC_DUMP_ROOT_FIXTURE
#endif // CbcRootFixtureDump_H
