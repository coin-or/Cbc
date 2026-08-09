/**
 * Fixture dumping for clique separation / clique extension experiments.
 *
 * @file CbcClqFixtureDump.hpp
 * @brief dump the exact inputs of a Bron-Kerbosch call, for offline replay
 *
 * Both users of Bron-Kerbosch inside CBC -- CglBKClique separating violated
 * clique inequalities, and CglCliqueStrengthening extending existing cliques --
 * only run deep inside a solve, after preprocessing and (for separation) after
 * the root LP. Experimenting on them therefore costs a full CBC run per
 * measurement, which is far too slow to iterate on the algorithm itself.
 *
 * This writes out everything such a call depends on, so a stand-alone tool can
 * reproduce it exactly:
 *
 *   <name>.<tag>.mps.gz    the problem as it stands at the call
 *   <name>.<tag>.cgraph    the conflict graph, serialized (not rebuilt: see below)
 *   <name>.<tag>.bas       the LP basis, so the optimum is one warm start away
 *   <name>.<tag>.sol       the LP solution, for cross-checking the warm start
 *   <name>.<tag>.meta      key/value provenance, read back by the benchmarks
 *
 * The `.meta` file carries the two fields a replay cannot infer: `lpOptimal` and
 * `extMethodRequested`. Clique extension runs twice in CbcSolver, and the two
 * passes are not interchangeable -- "before" runs on the original problem with no
 * LP solved and asks for extension method 2, "after" runs on the preprocessed one
 * with an optimal LP and asks for method 4. Method 4 needs reduced costs and
 * degrades to 2 without them, so a replay that defaulted to 4 everywhere would
 * measure neither pass faithfully: it would hit that degradation on the "before"
 * fixture, and there is no way to tell from the files alone which method was
 * wanted. Both facts are therefore recorded rather than guessed.
 *
 * The graph is saved rather than rebuilt on the other side because a rebuild
 * does not reproduce it. CBC's graph is built before CglCliqueStrengthening
 * deletes and adds rows, and OsiSolverInterface::checkCGraph() only rebuilds
 * when the *column* count changes -- so the graph a separator sees is routinely
 * stale with respect to the matrix it is separating on. That is deliberate, and
 * it is the graph whose behaviour we want to study, so it is the graph we store.
 * Measured: on trdtaunimep a rebuild gives 190900 direct conflicts against the
 * 191028 CBC actually used, and on co-100 3531734 against 3531732.
 *
 * Entirely behind CBC_DUMP_CLQSEP_FIXTURE and off by default: build with
 * -DCBC_DUMP_CLQSEP_FIXTURE to generate fixtures. Header-only static functions,
 * so no Makefile.am/Makefile.in changes are needed to carry it.
 *
 * Environment:
 *   CBC_CLQSEP_FIXTURE_DIR   output directory
 *                            (default ~/instances/miplib/2017+spp/clqsepFixtures)
 *   CBC_CLQSEP_FIXTURE_NAME  instance base name, when the MPS problem name is
 *                            unhelpful (drivers normally set this)
 **/

#ifndef CbcClqFixtureDump_H
#define CbcClqFixtureDump_H

#ifdef CBC_DUMP_CLQSEP_FIXTURE

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

#include "CoinStaticConflictGraph.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

/// mkdir -p, so a driver need not pre-create the tree.
static void cbcClqFixtureMkdirP(const std::string &path)
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
static std::string cbcClqFixtureDir()
{
  const char *env = getenv("CBC_CLQSEP_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/miplib/2017+spp/clqsepFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback,
/// but drivers should set CBC_CLQSEP_FIXTURE_NAME: several instances share a
/// problem name, and some carry none at all.
static std::string cbcClqFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_CLQSEP_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

/**
 * Write the problem, preserving the column *indices* the conflict graph refers to.
 *
 * CoinMpsIO::writeMps only emits a column that has matrix entries or a nonzero
 * objective ("only put out if elements or objective value"), so a column that is
 * empty in the matrix silently disappears. That is fatal here for two reasons,
 * the second worse than the first: the reloaded model has fewer columns than the
 * graph has nodes (n vs 2n bookkeeping breaks, and CglBKClique aborts outright),
 * and every column after the first dropped one shifts down an index, so the
 * graph's node numbering would silently point at the wrong variables.
 *
 * Measured on the original-problem captures: gmu-35-40 loses 363 of 1205 columns
 * and gmut-77-40 11198 of 24338, in both cases starting from the very first
 * column written -- so the shift is real, not hypothetical. Preprocessed captures
 * are unaffected, since preprocessing removes empty columns itself.
 *
 * The fix is one appended row holding a coefficient of 1 for every empty column,
 * which gives each of them a stored element and so gets them written at their
 * original index. It cannot be done with an explicit structural zero, which is
 * the obvious first idea and does not work: writeMps skips zero-valued elements
 * when emitting cards (`if (value && !stringRow[jRow])`), so a padded column
 * enters the column block and still writes nothing.
 *
 * Bounds on the pad row are the interval the sum of those columns can reach, so
 * the row is redundant and the LP is unchanged even for a consumer that never
 * removes it. When a column's bound is not finite that interval cannot be
 * computed, and 1e29 stands in -- just under CoinMpsIO's 1e30 infinity, so the
 * row still writes as a bounded 'L' rather than as a free 'N' row that the reader
 * would drop, taking the padding with it.
 *
 * Consumers should delete the row: `.meta` records the captured row count, so a
 * model that loads with exactly one row more than that ends in the pad row.
 * Round-trip checked on gmu-35-40 -- write, read back, delete the pad row, and
 * every column bound, objective coefficient, integrality flag, row bound and
 * matrix element matches the captured model, with the same 363 columns empty
 * again.
 *
 * Done on a clone, so the solver Cbc is in the middle of using is never touched.
 */
static int cbcClqFixtureWriteMps(const OsiSolverInterface *si, const std::string &path,
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

  if (empties.empty()) {
    // formatType 2 = IEEE hex, so every coefficient round-trips bit-for-bit;
    // CoinMpsIO adds free format itself when a name exceeds 8 characters.
    return si->writeMpsNative(path.c_str(), NULL, NULL, 2, 1);
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

  const int rc = clone->writeMpsNative(path.c_str(), NULL, NULL, 2, 1);
  delete clone;
  return rc;
}

/**
 * Write the basis. OsiClp's writeBasisNative() emits FREEIEEE, which round-trips
 * doubles exactly, and matches the .bas files of the existing
 * preProcessedInstances fixture set so the two remain interchangeable. The base
 * class implementation is a no-op that still returns success, so success is
 * confirmed by the file existing rather than by the return code; failing that,
 * fall back to a plain status dump good enough to rebuild a CoinWarmStartBasis.
 */
static bool cbcClqFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
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

/// Write the LP solution in the same shape as preProcessedInstances/*.sol:
/// an "=obj=" line, then one line per nonzero as "<index> <name> <value>".
static bool cbcClqFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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
 * Write the provenance file. `lpOptimal` and `extMethodRequested` are the
 * load-bearing fields: see the header comment. The rest is there so a fixture can
 * be read on its own, without the driver's CSV, and so a mismatch between graph
 * and matrix is visible without loading either.
 *
 * \param extMethodRequested the extension method the captured call site asked
 *        for, or -1 at a site where extension is not what runs (separation).
 */
static bool cbcClqFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  const CoinStaticConflictGraph *cgraph, int extMethodRequested, int paddedColumns,
  const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  fprintf(fp, "tag %s\n", tag);
  fprintf(fp, "rows %d\n", si->getNumRows());
  fprintf(fp, "cols %d\n", si->getNumCols());
  fprintf(fp, "elements %d\n", si->getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si->isProvenOptimal());
  fprintf(fp, "extMethodRequested %d\n", extMethodRequested);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index; see
  // cbcClqFixtureWriteMps. `rows` above is the captured count, so a consumer that
  // loads rows+1 should delete the last row to recover the captured model.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());
  fprintf(fp, "cgNodes %lu\n", (unsigned long)cgraph->size());
  fprintf(fp, "cgDirectConflicts %lu\n", (unsigned long)cgraph->nTotalDirectConflicts());
  fprintf(fp, "cgCliques %lu\n", (unsigned long)cgraph->nCliques());
  fprintf(fp, "cgCliqueElements %lu\n", (unsigned long)cgraph->nTotalCliqueElements());
  fprintf(fp, "cgDensity %.10f\n", cgraph->density());
  fprintf(fp, "cgMinDegree %lu\n", (unsigned long)cgraph->minDegree());
  fprintf(fp, "cgMaxDegree %lu\n", (unsigned long)cgraph->maxDegree());
  fprintf(fp, "cgUpdatedBounds %lu\n", (unsigned long)cgraph->updatedBounds().size());
  fprintf(fp, "cgInfeasibleImplications %lu\n",
    (unsigned long)cgraph->infeasibleImplications().size());
  fclose(fp);
  return true;
}

/**
 * Dump one fixture. Returns true if files were written.
 *
 * Instances whose graph holds no conflicts are skipped: neither BK consumer has
 * anything to do on them, so a fixture would only be noise. A large share of any
 * instance set falls in this bucket, which is why the skip is reported rather
 * than silent.
 *
 * \param si  solver whose matrix, LP solution and conflict graph are captured
 * \param tag distinguishes capture sites (e.g. "sep", "clqstr-before")
 * \param extMethodRequested extension method this site asked for; -1 where
 *        extension is not what runs, as at the separation site
 */
static bool cbcDumpClqFixture(OsiSolverInterface *si, const char *tag,
  int extMethodRequested = -1)
{
  if (!si)
    return false;

  const CoinStaticConflictGraph *cgraph = si->getCGraph();
  const std::string name = cbcClqFixtureName(si);

  if (!cgraph) {
    printf("[clqfixture] %s.%s: SKIP (no conflict graph)\n", name.c_str(), tag);
    fflush(stdout);
    return false;
  }
  if (cgraph->nTotalDirectConflicts() == 0 && cgraph->nTotalCliqueElements() == 0) {
    printf("[clqfixture] %s.%s: SKIP (graph has no conflicts; size=%lu)\n",
      name.c_str(), tag, (unsigned long)cgraph->size());
    fflush(stdout);
    return false;
  }

  const std::string dir = cbcClqFixtureDir();
  cbcClqFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcClqFixtureWriteMps(si, stem + ".mps", paddedColumns);
  const int cgRc = cgraph->save((stem + ".cgraph").c_str());
  const bool metaOk = cbcClqFixtureWriteMeta(si, tag, cgraph, extMethodRequested,
    paddedColumns, stem + ".meta");
  // Only when there is an LP to capture. In "before" mode there is not, and
  // writing an empty basis plus an obj=0 solution would look like a broken
  // fixture rather than a faithful one.
  const bool haveLp = si->isProvenOptimal();
  const bool basOk = haveLp ? cbcClqFixtureWriteBasis(si, stem + ".bas") : true;
  const bool solOk = haveLp ? cbcClqFixtureWriteSol(si, stem + ".sol") : true;

  printf("[clqfixture] %s.%s: %s rows=%d cols=%d padded=%d cgSize=%lu "
         "directConf=%lu cliques=%lu clqElems=%lu density=%.8f lpOptimal=%d "
         "obj=%.15g (mps=%d cgraph=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && cgRc == 0 && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), paddedColumns,
    (unsigned long)cgraph->size(),
    (unsigned long)cgraph->nTotalDirectConflicts(),
    (unsigned long)cgraph->nCliques(),
    (unsigned long)cgraph->nTotalCliqueElements(),
    cgraph->density(), (int)haveLp, si->getObjValue(),
    mpsRc, cgRc, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  return mpsRc == 0 && cgRc == 0;
}

#endif // CBC_DUMP_CLQSEP_FIXTURE
#endif // CbcClqFixtureDump_H
