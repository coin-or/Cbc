/**
 * Fixture dumping for CglProbing experiments.
 *
 * @file CbcProbingFixtureDump.hpp
 * @brief dump the exact inputs of a CglProbing call, for offline replay
 *
 * Same motivation as CbcClqFixtureDump.hpp and CbcZeroHalfFixtureDump.hpp:
 * CglProbing only runs deep inside a solve, after preprocessing and the root
 * LP, so experimenting on it costs a full CBC run per measurement. This
 * captures everything the call depends on:
 *
 *   <name>.probe.mps.gz   the problem as it stands at the call
 *   <name>.probe.bas      the LP basis, so the optimum is one warm start away
 *   <name>.probe.sol      the LP solution, for cross-checking the warm start
 *   <name>.probe.ctype    which columns are integer (MPS cannot say: see below)
 *   <name>.probe.meta     key/value provenance, read back by the benchmark
 *
 * WHY THE PAYLOAD IS AGAIN JUST THE GENERIC FIVE. CglProbing *does* carry
 * cached auxiliary structure that a snapshot would have to be serialized --
 * rowCopy_/columnCopy_/rowLower_/rowUpper_ from snapshot(), and the whole
 * clique block (numberCliques_, cliqueStart_, cliqueEntry_, oneFixStart_,
 * zeroFixStart_, endFixStart_, whichClique_, cliqueRow_, cliqueRowStart_) from
 * createCliques(). Neither is populated on CBC's path, which was confirmed by
 * reading rather than assumed:
 *
 *   - `grep -rn 'snapshot(' Cbc/src` finds nothing at all, so rowCopy_ is NULL
 *     and gutsOfGenerateCuts takes its "create from current" branch, building
 *     rowCopy from si.getMatrixByRow() on every call.
 *   - `createCliques` appears once in Cbc/src, at CbcModel.cpp:4827, and on a
 *     CglKnapsackCover -- never on a CglProbing. So numberCliques_ == 0 and the
 *     dispatch in gutsOfGenerateCuts goes to probe(), not probeCliques().
 *
 * So reconstructing from the .mps is exactly what CBC itself does, and there is
 * no analogue of the conflict-graph staleness trap. A bench that wants to
 * measure the snapshot or clique paths must call snapshot()/createCliques()
 * itself -- and should, because those are then *its* configuration, not
 * something the fixture silently decided.
 *
 * WHY NOT REUSE THE ZeroHalf OR CLIQUE FIXTURE POPULATION. The five *files* are
 * generator-agnostic and interchangeable; the *set of instances* is not. The
 * clique dump skips any instance whose conflict graph holds no conflicts; the
 * ZeroHalf dump skips any instance with no row of all-integer columns and
 * integral coefficients. Probing's precondition is neither: it needs an integer
 * column with a movable bound (`colUpper[i] - colLower[i] > 1.0e-8`), and cares
 * nothing about integrality of coefficients or about binary conflicts. A pure
 * general-integer model with fractional coefficients and no conflict graph is a
 * perfectly good, and often slow, Probing case. This dump gates on Probing's own
 * precondition -- see cbcProbeFixtureCount below.
 *
 * The `.ctype` sidecar exists because MPS cannot express "fixed and integer": a
 * column with lb == ub takes writeMps's " FX " branch, which has no integer
 * form, and reads back continuous. For Probing that is central rather than
 * cosmetic in two ways at once: `intVar[i]` decides whether a column is a
 * probing candidate at all, and it also decides whether a tightened bound
 * becomes an OsiColCut. A model whose fixed integers came back continuous would
 * both probe a different column set and emit different column cuts.
 *
 * WHAT THE .meta RECORDS ABOUT THE CALL, AND WHY EACH FIELD IS LOAD-BEARING.
 * Probing is the one generator CBC calls through a different entry point, so a
 * bench that just calls generateCuts() with a default CglTreeInfo measures
 * something CBC never runs. From CbcCutGenerator.cpp:340-376:
 *
 *   - CBC calls `generateCutsAndModify(solver, cs, info2)` where info2 is the
 *     model's CglTreeProbingInfo, not `generateCuts`. generateCutsAndModify
 *     leaves colLower_/colUpper_ behind for tightLower()/tightUpper().
 *   - `info2->options` is assigned outright as `globalCutsAtRoot() ? 8 : 0` --
 *     it is NOT or-ed with the bits the non-Probing branch builds up. So at the
 *     root, bits 64 (preprocessing), 2048 (playing around) and 16384 (just
 *     column cuts) are all clear. Consequences worth knowing before optimizing:
 *     probeFast() is unreachable from CBC's root call, and the FIXED_BOTH_WAYS
 *     blocks (gated on options&2048) never fire there either.
 *   - `info2->pass` is `model_->getCurrentPassNumber() - 1`, so it is 0 on the
 *     first root pass. That single fact is why pass 1 is the expensive one:
 *     probe() computes `justFix = (!inTree && !pass) ? -1 : 0` and then
 *     `if (justFix < 0) maxProbe = numberThisTime_`, which *overrides* the
 *     maxProbe==123 "be a bit intelligent" subsetting. Every candidate column
 *     is probed on pass 0; only later passes keep every 4th.
 *   - `info2->strengthenRow` is NULL on CBC's path -- the only setters in the
 *     tree are CglPreProcess.cpp:6749 and :7926. So row cuts land in `cs`
 *     rather than replacing rows, and the `needEffectiveness` threshold in
 *     probe() is 1.0e-3, not 1.0e-8.
 *   - `setMaxSeconds()` is applied only in CbcCutGenerator's *non*-Probing
 *     branch, so maxSeconds_ is 0 and probingDeadline is 0: the root Probing
 *     call has no internal time limit at all. `maxSeconds` in the .meta is the
 *     dumping run's budget, provenance only -- not a limit the call honoured.
 *
 * The knob values are recorded too, so a bench reproduces CBC's configuration
 * (CbcSolverCutSetup.cpp:90-127) instead of CglProbing's constructor defaults,
 * which differ on nearly every field.
 *
 * Entirely behind CBC_DUMP_PROBING_FIXTURE and off by default. Header-only
 * static functions, so no Makefile.am/Makefile.in changes are needed.
 *
 * Environment:
 *   CBC_PROBE_FIXTURE_DIR    output directory
 *                            (default ~/instances/miplib/2017+spp/probingFixtures)
 *   CBC_PROBE_FIXTURE_NAME   instance base name, when the MPS problem name is
 *                            unhelpful (drivers normally set this)
 **/

#ifndef CbcProbingFixtureDump_H
#define CbcProbingFixtureDump_H

#ifdef CBC_DUMP_PROBING_FIXTURE

#include <cmath>
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

#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

/// mkdir -p, so a driver need not pre-create the tree.
static void cbcProbeFixtureMkdirP(const std::string &path)
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
static std::string cbcProbeFixtureDir()
{
  const char *env = getenv("CBC_PROBE_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/miplib/2017+spp/probingFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback,
/// but drivers should set CBC_PROBE_FIXTURE_NAME: several instances share a
/// problem name, and some carry none at all.
static std::string cbcProbeFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_PROBE_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

/**
 * Write the MPS. formatType 2 is IEEE hex, which round-trips every double
 * exactly; plain MPS loses bits and the replayed LP then lands on a different
 * vertex -- which for Probing changes `lookedAt_`'s sort order and so the
 * probing order itself. numberAcross 1 for the same reason.
 *
 * Empty columns are the one wrinkle: writeMps drops a column with no matrix
 * entry, which shifts every later column index and so invalidates the .ctype
 * sidecar, the .bas and the .sol. A single redundant final row covering those
 * columns keeps them at their original index; `paddedColumns` in the .meta says
 * how many there were, and the loader deletes the extra row.
 */
static int cbcProbeFixtureWriteMps(const OsiSolverInterface *si,
  const std::string &path, int &paddedColumns)
{
  OsiSolverInterface *clone = si->clone();
  paddedColumns = 0;

  const CoinPackedMatrix *byCol = clone->getMatrixByCol();
  std::vector< int > empties;
  if (byCol) {
    const int *len = byCol->getVectorLengths();
    for (int j = 0; j < clone->getNumCols(); ++j) {
      if (len[j] == 0)
        empties.push_back(j);
    }
  }
  paddedColumns = (int)empties.size();

  std::vector< std::string > rowStore, colStore;
  std::vector< const char * > rowPtrs, colPtrs;

  if (!empties.empty()) {
    // A row that cannot cut anything off: its bounds are the interval the row
    // activity can take, so it is redundant by construction. Note this row is
    // free-ish rather than truly free, deliberately: gutsOfGenerateCuts deletes
    // rows with both bounds infinite, so a genuinely free pad row would vanish
    // and take the column indices with it.
    double padLb = 0.0, padUb = 0.0;
    const double *cl = clone->getColLower();
    const double *cu = clone->getColUpper();
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
  }

  // Names are regenerated so that the .bas written next matches them. A .bas
  // carrying names the .mps does not have is skipped by readBasis *and returns
  // success*, leaving the replay to start cold -- so this is load-bearing.
  const int nr = clone->getNumRows(), nc = clone->getNumCols();
  rowStore.reserve(nr);
  colStore.reserve(nc);
  for (int i = 0; i < nr; ++i)
    rowStore.push_back(clone->getRowName(i));
  for (int j = 0; j < nc; ++j)
    colStore.push_back(clone->getColName(j));
  rowPtrs.reserve(nr);
  colPtrs.reserve(nc);
  for (int i = 0; i < nr; ++i)
    rowPtrs.push_back(rowStore[i].c_str());
  for (int j = 0; j < nc; ++j)
    colPtrs.push_back(colStore[j].c_str());

  const int rc = clone->writeMpsNative(path.c_str(),
    rowPtrs.empty() ? NULL : &rowPtrs[0],
    colPtrs.empty() ? NULL : &colPtrs[0], 2, 1);
  delete clone;
  return rc;
}

/**
 * Write the basis. OsiClp's writeBasisNative() emits FREEIEEE, which round-trips
 * doubles exactly and matches the .bas files of preProcessedInstances so the two
 * sets stay interchangeable. The base-class implementation is a no-op that still
 * returns success, so success is confirmed by the file existing.
 */
static bool cbcProbeFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
{
  si->writeBasisNative(path.c_str());
  FILE *probe = fopen(path.c_str(), "r");
  if (probe) {
    fclose(probe);
    return true;
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
static bool cbcProbeFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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
 * `isContinuous()` is the source of truth rather than `getColType()`, which
 * derives Binary/GeneralInteger from the current bounds.
 */
static bool cbcProbeFixtureWriteColTypes(const OsiSolverInterface *si,
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
 * Count what CglProbing would actually work on, mirroring its own tests.
 *
 * `probeCandidates` is the acceptance test that decides whether the generator
 * can fire at all: gutsOfGenerateCuts fills lookedAt_ with every i where
 * `intVar[i] && colUpper[i] - colLower[i] > 1.0e-8`, and probe() then loops
 * `iLook` over exactly that. Zero candidates means the probing loop body never
 * executes -- nothing to measure. This is the skip gate.
 *
 * `probeBinaries` splits out the binary subset because the two kinds of
 * candidate cost differently: a binary is probed both ways from a 0/1 split,
 * while a general integer probes around floor/ceil of its LP value.
 *
 * `probeRows` / `probeRowNonzeros` reproduce the row filter, which is where the
 * propagation cost actually comes from, including its non-obvious escape hatch:
 * rows longer than maxElements or with both bounds infinite are deleted, BUT if
 * the deleted nonzeros are under a tenth of the matrix then `maxElements` is
 * reset to nCols and only the free rows go. So a model with a few dense rows
 * keeps them all, and the counts differ sharply from a naive length filter.
 * maxElementsRoot is used, matching a root call.
 */
static void cbcProbeFixtureCount(const OsiSolverInterface *si, int maxElementsRoot,
  int &probeCandidates, int &probeBinaries, int &probeRows, int &probeRowNonzeros)
{
  probeCandidates = probeBinaries = probeRows = probeRowNonzeros = 0;

  const int numberColumns = si->getNumCols();
  const int numberRows = si->getNumRows();
  const double *columnLower = si->getColLower();
  const double *columnUpper = si->getColUpper();
  const double *rowLower = si->getRowLower();
  const double *rowUpper = si->getRowUpper();
  const char *intVar = si->getColType(true);

  for (int j = 0; j < numberColumns; ++j) {
    if (!intVar[j] || si->isContinuous(j))
      continue;
    if (columnUpper[j] - columnLower[j] <= 1.0e-8)
      continue;
    ++probeCandidates;
    if (columnLower[j] == 0.0 && columnUpper[j] == 1.0)
      ++probeBinaries;
  }

  const CoinPackedMatrix *byRow = si->getMatrixByRow();
  if (!byRow)
    return;
  const int *rowLength = byRow->getVectorLengths();
  const CoinBigIndex nElements = byRow->getNumElements();

  // "keep all if only a few dense": gutsOfGenerateCuts' own rule.
  int maxElements = maxElementsRoot;
  CoinBigIndex nTotalOut = 0;
  for (int i = 0; i < numberRows; ++i) {
    if (rowLength[i] > maxElements
      || (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20))
      nTotalOut += rowLength[i];
  }
  if (nTotalOut * 10 < nElements)
    maxElements = numberColumns;

  for (int i = 0; i < numberRows; ++i) {
    if (rowLength[i] > maxElements
      || (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20))
      continue;
    ++probeRows;
    probeRowNonzeros += rowLength[i];
  }
}

/**
 * Write the provenance file.
 *
 * Beyond the usual model shape, this records the *call shape* -- see the header
 * comment for why each of pass/inTree/options/strengthenRow matters, and why a
 * bench that guesses them measures a call CBC never makes. `cutoff` is recorded
 * because usingObjective_ > 0 makes gutsOfGenerateCuts append the objective as
 * an extra row when the cutoff is finite, changing nRows and the row filter.
 */
static bool cbcProbeFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  int paddedColumns, int integerColumns, int probeCandidates, int probeBinaries,
  int probeRows, int probeRowNonzeros, int pass, int options, int inTree,
  int formulationRows, int strengthenRow, double cutoff, double maxSeconds,
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
  // The dumping run's budget, NOT a limit the Probing call honoured: CBC never
  // calls setMaxSeconds on a CglProbing, so probingDeadline is 0 there.
  fprintf(fp, "maxSeconds %.6f\n", maxSeconds);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index. `rows` above
  // is the captured count, so a consumer loading rows+1 deletes the last row.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());
  // What CglProbing would look at: see cbcProbeFixtureCount.
  fprintf(fp, "probeCandidates %d\n", probeCandidates);
  fprintf(fp, "probeBinaries %d\n", probeBinaries);
  fprintf(fp, "probeRows %d\n", probeRows);
  fprintf(fp, "probeRowNonzeros %d\n", probeRowNonzeros);
  // The call shape CBC used, so a bench can reproduce it rather than guess.
  fprintf(fp, "infoPass %d\n", pass);
  fprintf(fp, "infoOptions %d\n", options);
  fprintf(fp, "infoInTree %d\n", inTree);
  fprintf(fp, "infoFormulationRows %d\n", formulationRows);
  fprintf(fp, "infoStrengthenRow %d\n", strengthenRow);
  fprintf(fp, "cutoff %.15g\n", cutoff);
  // CBC's Probing configuration (CbcSolverCutSetup.cpp:90-127), which differs
  // from CglProbing's constructor defaults on nearly every field. Recorded so
  // the bench's defaults can be CBC's, not the library's.
  fprintf(fp, "cbcMode 1\n");
  fprintf(fp, "cbcUsingObjective 1\n");
  fprintf(fp, "cbcMaxPassRoot 1\n");
  fprintf(fp, "cbcMaxProbeRoot 123\n");
  fprintf(fp, "cbcMaxLookRoot 20\n");
  fprintf(fp, "cbcMaxElementsRoot 300\n");
  fprintf(fp, "cbcRowCuts -3\n");
  fclose(fp);
  return true;
}

/**
 * Dump one fixture. Returns true if files were written.
 *
 * The skip criterion is Probing's own: with no integer column of movable bounds
 * there is no probing candidate, probe()'s `iLook` loop body never runs, and
 * there is nothing to measure. Deliberately *not* gated on a conflict graph or
 * on integral coefficients -- see the header comment on why inheriting either
 * population would bias this.
 *
 * A fixture without a proven-optimal LP is still written, but `lpOptimal 0` in
 * the .meta warns the consumer that the .bas/.sol are absent and the replay must
 * solve from cold.
 */
static bool cbcDumpProbingFixture(OsiSolverInterface *si, const char *tag = "probe",
  int pass = 0, int options = 0, int inTree = 0, int formulationRows = -1,
  int strengthenRow = 0, double maxSeconds = 0.0)
{
  if (!si)
    return false;

  const std::string name = cbcProbeFixtureName(si);

  if (si->getNumCols() == 0 || si->getNumRows() == 0) {
    printf("[probefixture] %s.%s: SKIP (empty model: %d rows, %d cols)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols());
    fflush(stdout);
    return false;
  }

  int probeCandidates = 0, probeBinaries = 0, probeRows = 0, probeRowNonzeros = 0;
  // 300 is CBC's maxElementsRoot; see cbcProbeFixtureWriteMeta.
  cbcProbeFixtureCount(si, 300, probeCandidates, probeBinaries, probeRows,
    probeRowNonzeros);
  if (probeCandidates == 0) {
    printf("[probefixture] %s.%s: SKIP (no probing candidate; %d rows, %d cols "
           "-- no integer column has colUpper-colLower > 1e-8, so probe()'s "
           "iLook loop is empty)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols());
    fflush(stdout);
    return false;
  }

  double cutoff = COIN_DBL_MAX;
  si->getDblParam(OsiDualObjectiveLimit, cutoff);

  const std::string dir = cbcProbeFixtureDir();
  cbcProbeFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcProbeFixtureWriteMps(si, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk = cbcProbeFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const bool metaOk = cbcProbeFixtureWriteMeta(si, tag, paddedColumns,
    integerColumns, probeCandidates, probeBinaries, probeRows, probeRowNonzeros,
    pass, options, inTree,
    formulationRows >= 0 ? formulationRows : si->getNumRows(),
    strengthenRow, cutoff, maxSeconds, stem + ".meta");
  const bool haveLp = si->isProvenOptimal();
  const bool basOk = haveLp ? cbcProbeFixtureWriteBasis(si, stem + ".bas") : true;
  const bool solOk = haveLp ? cbcProbeFixtureWriteSol(si, stem + ".sol") : true;

  printf("[probefixture] %s.%s: %s rows=%d cols=%d int=%d padded=%d "
         "cand=%d bin=%d probeRows=%d probeNz=%d lpOptimal=%d obj=%.15g "
         "pass=%d options=%d (mps=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, paddedColumns,
    probeCandidates, probeBinaries, probeRows, probeRowNonzeros,
    (int)haveLp, si->getObjValue(), pass, options,
    mpsRc, (int)ctOk, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  return mpsRc == 0 && ctOk && metaOk;
}

#endif // CBC_DUMP_PROBING_FIXTURE
#endif // CbcProbingFixtureDump_H
