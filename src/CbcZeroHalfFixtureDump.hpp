/**
 * Fixture dumping for 0-1/2 (ZeroHalf) cut separation experiments.
 *
 * @file CbcZeroHalfFixtureDump.hpp
 * @brief dump the exact inputs of a CglZeroHalf separation call, for offline replay
 *
 * Same motivation as CbcClqFixtureDump.hpp: CglZeroHalf only runs deep inside a
 * solve, after preprocessing and the root LP, so experimenting on it costs a full
 * CBC run per measurement. This captures everything the call depends on:
 *
 *   <name>.zh.mps.gz   the problem as it stands at the call
 *   <name>.zh.bas      the LP basis, so the optimum is one warm start away
 *   <name>.zh.sol      the LP solution, for cross-checking the warm start
 *   <name>.zh.ctype    which columns are integer (MPS cannot say: see below)
 *   <name>.zh.meta     key/value provenance, read back by the benchmark
 *
 * The payload is deliberately the *generator-agnostic* five and nothing else.
 * CglZeroHalf has no cached auxiliary structure to serialize: refreshSolver()
 * derives mr_/mc_/mnz_/mtbeg_/mtcnt_/mtind_/mtval_/vlb_/vub_/mrhs_/msense_ from
 * the solver's row-major matrix, column bounds and column types every time it is
 * called, and generateCuts() re-reads the bounds on top of that. So there is no
 * analogue of the conflict-graph staleness trap here -- reconstructing on load is
 * not merely acceptable, it is exactly what CBC itself does. Verified by reading
 * refreshSolver(): every field it fills comes from the solver, none survives from
 * a previous model.
 *
 * WHY THIS EXISTS ALONGSIDE THE CLIQUE FIXTURES. The 237 `.sep` fixtures were
 * selected by a clique criterion -- cbcDumpClqFixture skips any instance whose
 * conflict graph holds no conflicts, and one whose graph does not span 2n nodes.
 * That is the wrong population for ZeroHalf, which keys off *integer
 * coefficients* and cares nothing about binary conflicts: an instance of pure
 * general integers with no conflict graph at all is a perfectly good, and
 * possibly hard, ZeroHalf case. Reusing the clique fixture set would therefore
 * silently exclude 121 of the 358 instances and bias every measurement. This
 * dump gates on ZeroHalf's own precondition instead (see below).
 *
 * The `.ctype` sidecar exists because MPS cannot express "fixed and integer":
 * a column with lb == ub takes writeMps's " FX " branch, which has no integer
 * form, and comes back continuous. For ZeroHalf that is not cosmetic but
 * central -- refreshSolver() rejects any row touching a continuous column
 * (CglZeroHalf.cpp: `if (vlb_[jColumn]==COIN_INT_MAX) good=false`), so every
 * integer column CBC had fixed by bound tightening would delete every row it
 * appears in, shrinking mr_/mnz_ and changing what is separated.
 *
 * Entirely behind CBC_DUMP_ZEROHALF_FIXTURE and off by default. Header-only
 * static functions, so no Makefile.am/Makefile.in changes are needed.
 *
 * Environment:
 *   CBC_ZH_FIXTURE_DIR    output directory
 *                         (default ~/instances/miplib/2017+spp/zerohalfFixtures)
 *   CBC_ZH_FIXTURE_NAME   instance base name, when the MPS problem name is
 *                         unhelpful (drivers normally set this)
 **/

#ifndef CbcZeroHalfFixtureDump_H
#define CbcZeroHalfFixtureDump_H

#ifdef CBC_DUMP_ZEROHALF_FIXTURE

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
static void cbcZhFixtureMkdirP(const std::string &path)
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
static std::string cbcZhFixtureDir()
{
  const char *env = getenv("CBC_ZH_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/miplib/2017+spp/zerohalfFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback,
/// but drivers should set CBC_ZH_FIXTURE_NAME: several instances share a problem
/// name, and some carry none at all.
static std::string cbcZhFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_ZH_FIXTURE_NAME");
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
 * vertex. numberAcross 1 for the same reason.
 *
 * Empty columns are the one wrinkle: writeMps drops a column with no matrix
 * entry, which shifts every later column index and so invalidates the .ctype
 * sidecar, the .bas and the .sol. A single redundant final row covering those
 * columns keeps them at their original index; `paddedColumns` in the .meta says
 * how many there were, and the loader deletes the extra row.
 */
static int cbcZhFixtureWriteMps(const OsiSolverInterface *si,
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
    // activity can take, so it is redundant by construction.
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
static bool cbcZhFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
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
static bool cbcZhFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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
static bool cbcZhFixtureWriteColTypes(const OsiSolverInterface *si,
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
 * Count the rows CglZeroHalf would actually keep, mirroring refreshSolver()'s
 * acceptance test exactly: every column in the row integer, every coefficient
 * integral and below COIN_INT_MAX, both row bounds integral, and |rhs| in range.
 * A ranged row counts twice, as refreshSolver() splits it into a 'G' and an 'L'.
 *
 * This is recorded in the .meta so a driver can tell "ZeroHalf has nothing to do
 * here" from "the fixture failed to capture", and can filter the set down to the
 * instances that actually exercise the generator without loading each model.
 */
static void cbcZhFixtureCountIlp(const OsiSolverInterface *si,
  int &zhRows, int &zhNonzeros, int &zhIntegerCols)
{
  zhRows = zhNonzeros = zhIntegerCols = 0;

  const CoinPackedMatrix *byRow = si->getMatrixByRow();
  if (!byRow)
    return;
  const int *column = byRow->getIndices();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *rowLength = byRow->getVectorLengths();
  const double *rowElements = byRow->getElements();
  const double *columnLower = si->getColLower();
  const double *columnUpper = si->getColUpper();
  const double *rowLower = si->getRowLower();
  const double *rowUpper = si->getRowUpper();
  const char *intVar = si->getColType(true);

  const int numberColumns = si->getNumCols();
  const int numberRows = si->getNumRows();

  // Reproduce refreshSolver's vlb_ sentinel: COIN_INT_MAX marks "continuous",
  // which is what makes a row unusable.
  std::vector< int > vlb(numberColumns);
  for (int j = 0; j < numberColumns; ++j) {
    if (intVar[j]) {
      double lo = columnLower[j];
      if (lo < -COIN_INT_MAX)
        lo = -COIN_INT_MAX;
      vlb[j] = (int)ceil(lo);
      ++zhIntegerCols;
    } else {
      vlb[j] = COIN_INT_MAX;
    }
  }

  for (int i = 0; i < numberRows; ++i) {
    const int n = rowLength[i];
    bool good = (n > 0);
    for (CoinBigIndex j = rowStart[i]; j < rowStart[i] + n; ++j) {
      if (vlb[column[j]] == COIN_INT_MAX) {
        good = false;
        break;
      }
      const double value = rowElements[j];
      if (fabs(value - floor(value + 0.5)) > 1.0e-15 || fabs(value) >= COIN_INT_MAX) {
        good = false;
        break;
      }
    }
    const double lo = rowLower[i], up = rowUpper[i];
    int iType = 1;
    double rhs = 1.0e20;
    if (lo > -1.0e20) {
      if (fabs(lo - floor(lo + 0.5)) > 1.0e-15)
        good = false;
      rhs = fabs(lo);
      if (up < 1.0e20) {
        rhs = std::max(fabs(lo), fabs(up));
        if (lo != up)
          iType = 2;
        if (fabs(up - floor(up + 0.5)) > 1.0e-12)
          good = false;
      }
    } else if (up < 1.0e20) {
      rhs = fabs(up);
      if (fabs(up - floor(up + 0.5)) > 1.0e-12)
        good = false;
    }
    if (good && rhs < COIN_INT_MAX) {
      zhRows += iType;
      zhNonzeros += iType * n;
    }
  }
}

/// Write the provenance file. `lpOptimal` and the `zh*` counts are the
/// load-bearing fields: see cbcZhFixtureCountIlp.
static bool cbcZhFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  int paddedColumns, int integerColumns, int zhRows, int zhNonzeros,
  int zhIntegerCols, double maxSeconds, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;
  fprintf(fp, "tag %s\n", tag);
  fprintf(fp, "rows %d\n", si->getNumRows());
  fprintf(fp, "cols %d\n", si->getNumCols());
  fprintf(fp, "elements %d\n", si->getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si->isProvenOptimal());
  fprintf(fp, "maxSeconds %.6f\n", maxSeconds);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index. `rows` above
  // is the captured count, so a consumer loading rows+1 deletes the last row.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());
  // What CglZeroHalf::refreshSolver() would build from this model.
  fprintf(fp, "zhRows %d\n", zhRows);
  fprintf(fp, "zhNonzeros %d\n", zhNonzeros);
  fprintf(fp, "zhIntegerCols %d\n", zhIntegerCols);
  fclose(fp);
  return true;
}

/**
 * Dump one fixture. Returns true if files were written.
 *
 * The skip criterion is ZeroHalf's own: if refreshSolver() would find no usable
 * row (zhNonzeros == 0) then generateCuts() returns immediately on the `mnz_`
 * test and there is nothing to measure. Deliberately *not* gated on anything to
 * do with a conflict graph -- see the header comment on why reusing the clique
 * fixture population would bias this.
 *
 * A fixture without a proven-optimal LP is still written, but `lpOptimal 0` in
 * the .meta warns the consumer that the .bas/.sol are absent and the replay must
 * solve from cold.
 */
static bool cbcDumpZeroHalfFixture(OsiSolverInterface *si, const char *tag = "zh",
  double maxSeconds = 0.0)
{
  if (!si)
    return false;

  const std::string name = cbcZhFixtureName(si);

  if (si->getNumCols() == 0 || si->getNumRows() == 0) {
    printf("[zhfixture] %s.%s: SKIP (empty model: %d rows, %d cols)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols());
    fflush(stdout);
    return false;
  }

  int zhRows = 0, zhNonzeros = 0, zhIntegerCols = 0;
  cbcZhFixtureCountIlp(si, zhRows, zhNonzeros, zhIntegerCols);
  if (zhNonzeros == 0) {
    printf("[zhfixture] %s.%s: SKIP (no ZeroHalf-usable row; %d rows, %d cols, "
           "%d integer cols -- generateCuts would return on mnz_==0)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols(), zhIntegerCols);
    fflush(stdout);
    return false;
  }

  const std::string dir = cbcZhFixtureDir();
  cbcZhFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcZhFixtureWriteMps(si, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk = cbcZhFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const bool metaOk = cbcZhFixtureWriteMeta(si, tag, paddedColumns, integerColumns,
    zhRows, zhNonzeros, zhIntegerCols, maxSeconds, stem + ".meta");
  const bool haveLp = si->isProvenOptimal();
  const bool basOk = haveLp ? cbcZhFixtureWriteBasis(si, stem + ".bas") : true;
  const bool solOk = haveLp ? cbcZhFixtureWriteSol(si, stem + ".sol") : true;

  printf("[zhfixture] %s.%s: %s rows=%d cols=%d int=%d padded=%d zhRows=%d "
         "zhNz=%d zhIntCols=%d lpOptimal=%d obj=%.15g "
         "(mps=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, paddedColumns,
    zhRows, zhNonzeros, zhIntegerCols, (int)haveLp, si->getObjValue(),
    mpsRc, (int)ctOk, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  return mpsRc == 0 && ctOk && metaOk;
}

#endif // CBC_DUMP_ZEROHALF_FIXTURE
#endif // CbcZeroHalfFixtureDump_H
