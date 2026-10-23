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
 *   <name>.<tag>.ctype     which columns are integer (MPS cannot say: see below)
 *   <name>.<tag>.meta      key/value provenance, read back by the benchmarks
 *
 * The `.ctype` sidecar exists because MPS cannot express "fixed and integer".
 * CoinMpsIO::writeMps carries integrality purely in the bound *type* -- BV, UI,
 * LI, or MI+UI -- and a column with lb == ub takes the " FX " branch
 * (CoinMpsIO.cpp:4725), which has no integer form; the reader in turn only marks
 * a column integer from an INTORG marker or one of those bound types, and
 * writeMps emits no markers at all. So every integer column CBC had *fixed* by
 * bound tightening comes back continuous. That is not cosmetic: both BK
 * consumers gate on integrality, CglCliqueStrengthening::detectCliqueRows most
 * visibly, since it rejects any row holding a non-binary column. Measured on
 * physiciansched3-3: 39994 of 72141 integer columns lost -- every one of them
 * fixed -- taking clique rows from 188486 down to 83611, and with them 1759 of
 * the extensions CBC found. The sidecar records CoinColumnType per column and is
 * re-applied on load, restoring the captured model exactly.
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
 * **Names must be passed explicitly, not left to CoinMpsIO.** Passing NULL for
 * `columnNames` makes it generate `C0000000`/`R0000000`, while
 * ClpSimplexOther::writeBasis prints the model's *real* names -- so the `.bas`
 * says `x(0,0)` and the MPS says `C0000000`. CoinMpsIO::readBasis resolves every
 * line by name lookup (`findHash`) and silently skips what it cannot find, so
 * every column missed and the reader was left with an all-slack basis. That reads
 * as success -- rc=1, no warning -- and then costs a full cold solve: measured on
 * decomp2, `basicCols=0` and 3390 resolve iterations against 0 with real names,
 * landing on a different vertex with 1123 columns moved. Since BK separates
 * cliques violated by the *current fractional point*, that silently changed which
 * cliques the replay found. The pre-existing `preProcessedInstances/*.bas` set has
 * the same defect, for the same reason.
 *
 * Done on a clone, so the solver Cbc is in the middle of using is never touched.
 */
/// Collect the model's row and column names into the `const char **` pair
/// writeMpsNative wants, so the MPS carries the same names the `.bas` will.
/// getRowName/getColName are non-const, hence the non-const solver argument.
static void cbcClqFixtureNames(OsiSolverInterface *si,
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

static int cbcClqFixtureWriteMps(OsiSolverInterface *si, const std::string &path,
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
    cbcClqFixtureNames(si, rowStore, colStore, rowPtrs, colPtrs);
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

  // Named after the pad row is added, so the extra row gets a name too; the
  // consumer deletes it anyway.
  cbcClqFixtureNames(clone, rowStore, colStore, rowPtrs, colPtrs);
  const int rc = clone->writeMpsNative(path.c_str(),
    rowPtrs.empty() ? NULL : &rowPtrs[0],
    colPtrs.empty() ? NULL : &colPtrs[0], 2, 1);
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
 * Write the column types, which the MPS cannot carry for fixed integers.
 *
 * Only the integer columns are listed, as "<index> <CoinColumnType>", since they
 * are usually the minority and the default on load is continuous. A header line
 * gives the column count so a loader can reject a sidecar written for a
 * different model. `isContinuous()` is the source of truth here rather than
 * `getColType()`, which derives Binary/GeneralInteger from the current bounds --
 * that derivation is fine to reproduce on load, but it is not what needs saving.
 */
static bool cbcClqFixtureWriteColTypes(const OsiSolverInterface *si,
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
 * Write the provenance file. `lpOptimal` and `extMethodRequested` are the
 * load-bearing fields: see the header comment. The rest is there so a fixture can
 * be read on its own, without the driver's CSV, and so a mismatch between graph
 * and matrix is visible without loading either.
 *
 * \param extMethodRequested the extension method the captured call site asked
 *        for, or -1 at a site where extension is not what runs (separation).
 * \param maxSeconds the wall-clock limit the captured site imposed, or 0 for
 *        none. Recorded because it changes the outcome and cannot be inferred:
 *        clique strengthening gets whatever is left of the model's budget, so a
 *        replay with no limit does more work than the capture did.
 */
static bool cbcClqFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  const CoinStaticConflictGraph *cgraph, int extMethodRequested, int paddedColumns,
  int integerColumns, double maxSeconds, const std::string &path)
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
  fprintf(fp, "maxSeconds %.6f\n", maxSeconds);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index; see
  // cbcClqFixtureWriteMps. `rows` above is the captured count, so a consumer that
  // loads rows+1 should delete the last row to recover the captured model.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  // How many columns the .ctype sidecar lists. A loader that restores fewer than
  // this has a stale or truncated sidecar, which would silently change which rows
  // count as cliques -- so it is worth being able to check.
  fprintf(fp, "integerColumns %d\n", integerColumns);
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
 * So are captures whose graph does not span the matrix's column space, which
 * cannot be replayed at all: the conflict graph is 2n nodes for n columns, and a
 * consumer that loaded a graph of the wrong size would be indexing arbitrary
 * variables. Normally the invariant holds by construction, since checkCGraph()
 * rebuilds whenever the column count changes -- but a model can also be *emptied*
 * after its graph was built, without that path being taken. Measured: wqueens-50,
 * -100, -200, -300 and h80x6320 all reach the "after" strengthening site with a
 * 0-row 0-column preprocessed model ("Processed model: 0 rows, 0 cols, 0 NZ" --
 * presolve finished them) while still carrying the original graph, so the dumped
 * fixture pairs a 5000-node graph with a 0-column matrix.
 *
 * \param si  solver whose matrix, LP solution and conflict graph are captured
 * \param tag distinguishes capture sites (e.g. "sep", "clqstr-before")
 * \param extMethodRequested extension method this site asked for; -1 where
 *        extension is not what runs, as at the separation site
 * \param maxSeconds wall-clock limit this site imposed; 0 for none
 */
static bool cbcDumpClqFixture(OsiSolverInterface *si, const char *tag,
  int extMethodRequested = -1, double maxSeconds = 0.0)
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
  if (cgraph->size() != (size_t)si->getNumCols() * 2) {
    printf("[clqfixture] %s.%s: SKIP (graph spans %lu nodes, model has %d cols so "
           "%d expected; not replayable)\n",
      name.c_str(), tag, (unsigned long)cgraph->size(), si->getNumCols(),
      si->getNumCols() * 2);
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
  int integerColumns = 0;
  const bool ctOk = cbcClqFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const bool metaOk = cbcClqFixtureWriteMeta(si, tag, cgraph, extMethodRequested,
    paddedColumns, integerColumns, maxSeconds, stem + ".meta");
  // Only when there is an LP to capture. In "before" mode there is not, and
  // writing an empty basis plus an obj=0 solution would look like a broken
  // fixture rather than a faithful one.
  const bool haveLp = si->isProvenOptimal();
  const bool basOk = haveLp ? cbcClqFixtureWriteBasis(si, stem + ".bas") : true;
  const bool solOk = haveLp ? cbcClqFixtureWriteSol(si, stem + ".sol") : true;

  printf("[clqfixture] %s.%s: %s rows=%d cols=%d int=%d padded=%d cgSize=%lu "
         "directConf=%lu cliques=%lu clqElems=%lu density=%.8f lpOptimal=%d "
         "obj=%.15g (mps=%d cgraph=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && cgRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, paddedColumns,
    (unsigned long)cgraph->size(),
    (unsigned long)cgraph->nTotalDirectConflicts(),
    (unsigned long)cgraph->nCliques(),
    (unsigned long)cgraph->nTotalCliqueElements(),
    cgraph->density(), (int)haveLp, si->getObjValue(),
    mpsRc, cgRc, (int)ctOk, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  return mpsRc == 0 && cgRc == 0 && ctOk;
}

#endif // CBC_DUMP_CLQSEP_FIXTURE
#endif // CbcClqFixtureDump_H
