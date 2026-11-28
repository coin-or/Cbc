/**
 * Fixture dumping for RINS (and other incumbent-based) heuristic experiments.
 *
 * @file CbcHeuristicRinsFixtureDump.hpp
 * @brief dump the current node's LP relaxation + the incumbent solution, for
 *        offline replay of RINS/VND-style "fix by agreement with the
 *        incumbent, then re-solve" heuristics.
 *
 * RINS (CbcHeuristicRINS::solution()) only becomes interesting to experiment
 * on once TWO things exist simultaneously: a fractional LP relaxation at some
 * node, and an incumbent integer-feasible solution to compare it against.
 * Reaching that state costs a real branch-and-bound run per experiment --
 * often much more than reaching the root LP alone, since an incumbent has to
 * be found first. This writes out everything needed to reconstruct that exact
 * state offline:
 *
 *   <name>.rins.mps.gz     the current node's matrix (preprocessed problem +
 *                          any cuts added so far), formatType 2 (IEEE hex).
 *                          Bounds reflect this node's branching decisions, so
 *                          a replay driver's fresh CbcModel will treat THIS
 *                          node's bounds as "original" -- see the caveat below.
 *   <name>.rins.bas        the optimal LP basis at this node (warm start,
 *                          0 iterations)
 *   <name>.rins.sol        the fractional LP solution at this node
 *   <name>.rins.ctype      integer/continuous per column (MPS cannot express a
 *                          fixed integer column)
 *   <name>.rins.incumbent  the incumbent (best-known) solution's value for
 *                          EVERY column, dense -- RINS needs the exact value,
 *                          including zeros, to decide agreement
 *   <name>.rins.meta       rows/cols/elements/nodeCount/depth/lpObjValue/
 *                          incumbentObjValue/numberIntegers, so a driver can
 *                          filter fixtures without loading them
 *
 * Caveat (documented, not fixed here): CbcModel derives each integer object's
 * "original bounds" from the solver's bounds at CONSTRUCTION time. The real
 * RINS call this fixture captures reads original bounds set at the ROOT
 * (unaffected by later branching), which affects `shallowDepth_` fixing modes
 * 1-3 (fix only at/away-from-original-lower-bound). A fresh CbcModel built
 * around this fixture instead treats the CAPTURED node's (branching-
 * tightened) bounds as original. For `shallowDepth_ == 0` (the CLI default,
 * "fix all agreeing") this makes no difference at all. For modes 1-3 it can
 * shift results: `mip-heur-replay`'s output notes this whenever a non-zero
 * shallow depth is requested, so it isn't silently misread as gospel.
 *
 * Entirely behind CBC_DUMP_RINS_FIXTURE and off by default: build with
 * -DCBC_DUMP_RINS_FIXTURE to generate fixtures. Header-only static functions,
 * so no Makefile.am/Makefile.in changes are needed to carry it -- same
 * approach as CbcRootFixtureDump.hpp (whose writer functions this
 * intentionally duplicates rather than shares, since that header's contents
 * are entirely gated behind its own -DCBC_DUMP_ROOT_FIXTURE macro).
 *
 * Environment:
 *   CBC_RINS_FIXTURE_DIR   output directory
 *                          (default ~/instances/mip-sanity-data/rinsFixtures)
 *   CBC_RINS_FIXTURE_NAME  instance base name, when the MPS problem name is
 *                          unhelpful (drivers normally set this)
 **/

#ifndef CbcHeuristicRinsFixtureDump_H
#define CbcHeuristicRinsFixtureDump_H

#ifdef CBC_DUMP_RINS_FIXTURE

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

#include "CbcModel.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

/// mkdir -p, so a driver need not pre-create the tree.
static void cbcRinsFixtureMkdirP(const std::string &path)
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
static std::string cbcRinsFixtureDir()
{
  const char *env = getenv("CBC_RINS_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/mip-sanity-data/rinsFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback,
/// but drivers should set CBC_RINS_FIXTURE_NAME.
static std::string cbcRinsFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_RINS_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

static void cbcRinsFixtureNames(OsiSolverInterface *si,
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
 * Write the problem, preserving column *indices* via the same empty-column
 * padding trick as CbcRootFixtureDump.hpp (see that file's header comment for
 * the measured rationale): CoinMpsIO::writeMps drops a column with no matrix
 * elements and a zero objective, which would otherwise shift every later
 * column index.
 */
static int cbcRinsFixtureWriteMps(OsiSolverInterface *si, const std::string &path,
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
    cbcRinsFixtureNames(si, rowStore, colStore, rowPtrs, colPtrs);
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

  cbcRinsFixtureNames(clone, rowStore, colStore, rowPtrs, colPtrs);
  const int rc = clone->writeMpsNative(path.c_str(),
    rowPtrs.empty() ? NULL : &rowPtrs[0],
    colPtrs.empty() ? NULL : &colPtrs[0], 2, 1);
  delete clone;
  return rc;
}

static bool cbcRinsFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
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

static bool cbcRinsFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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

static bool cbcRinsFixtureWriteColTypes(const OsiSolverInterface *si,
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

/// Write the incumbent's value for EVERY column, dense -- RINS needs exact
/// values, including zeros, to test agreement with the LP relaxation.
static bool cbcRinsFixtureWriteIncumbent(const double *bestSolution, int numberColumns,
  double incumbentObj, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  fprintf(fp, "cols %d\n", numberColumns);
  fprintf(fp, "obj %.15g\n", incumbentObj);
  for (int j = 0; j < numberColumns; ++j)
    fprintf(fp, "%d %.17g\n", j, bestSolution[j]);
  fclose(fp);
  return true;
}

static bool cbcRinsFixtureWriteMeta(const OsiSolverInterface *si, CbcModel *model,
  int paddedColumns, int integerColumns, double incumbentObj, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  fprintf(fp, "tag rins\n");
  fprintf(fp, "rows %d\n", si->getNumRows());
  fprintf(fp, "cols %d\n", si->getNumCols());
  fprintf(fp, "elements %d\n", si->getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si->isProvenOptimal());
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "lpObjValue %.15g\n", si->getObjValue());
  fprintf(fp, "incumbentObjValue %.15g\n", incumbentObj);
  fprintf(fp, "nodeCount %lld\n", (long long)model->getNodeCount());
  fprintf(fp, "numberIntegers %d\n", model->numberIntegers());
  fclose(fp);
  return true;
}

/**
 * Dump one RINS fixture from the exact call site inside
 * CbcHeuristicRINS::solution(), right after the node/frequency gate passes
 * and before any variable is fixed. Returns true if files were written.
 *
 * \param model         the CbcModel driving the search (for node count /
 *                      numberIntegers)
 * \param solver        model_->solver(), holding the current node's LP
 *                      relaxation (already optimal at this call site)
 * \param bestSolution  the incumbent, one value per column
 */
static bool cbcDumpRinsFixture(CbcModel *model, OsiSolverInterface *solver,
  const double *bestSolution)
{
  if (!model || !solver || !bestSolution)
    return false;

  const std::string name = cbcRinsFixtureName(solver);
  const std::string dir = cbcRinsFixtureDir();
  cbcRinsFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + ".rins";

  int paddedColumns = 0;
  const int mpsRc = cbcRinsFixtureWriteMps(solver, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk = cbcRinsFixtureWriteColTypes(solver, stem + ".ctype", integerColumns);
  const bool basOk = cbcRinsFixtureWriteBasis(solver, stem + ".bas");
  const bool solOk = cbcRinsFixtureWriteSol(solver, stem + ".sol");

  // The incumbent's objective: recompute from bestSolution rather than trust
  // model->getObjValue(), which can reflect a LATER, better incumbent than
  // the one RINS is actually comparing against at this exact call. Raw dot
  // product plus offset, same convention getObjValue() itself reports in.
  double incumbentObj = 0.0;
  {
    const double *obj = solver->getObjCoefficients();
    const int n = solver->getNumCols();
    for (int j = 0; j < n; ++j)
      incumbentObj += obj[j] * bestSolution[j];
  }
  const bool incOk = cbcRinsFixtureWriteIncumbent(bestSolution, solver->getNumCols(),
    incumbentObj, stem + ".incumbent");
  const bool metaOk = cbcRinsFixtureWriteMeta(solver, model, paddedColumns, integerColumns,
    incumbentObj, stem + ".meta");

  printf("[rinsfixture] %s: %s rows=%d cols=%d int=%d padded=%d node=%lld lpObj=%.15g "
         "incumbentObj=%.15g (mps=%d ctype=%d bas=%d sol=%d inc=%d meta=%d)\n",
    name.c_str(),
    (mpsRc == 0 && ctOk && basOk && solOk && incOk && metaOk) ? "DUMPED" : "PARTIAL",
    solver->getNumRows(), solver->getNumCols(), integerColumns, paddedColumns,
    (long long)model->getNodeCount(), solver->getObjValue(), incumbentObj,
    mpsRc, (int)ctOk, (int)basOk, (int)solOk, (int)incOk, (int)metaOk);
  fflush(stdout);

  return mpsRc == 0 && ctOk && basOk && solOk && incOk && metaOk;
}

#endif // CBC_DUMP_RINS_FIXTURE
#endif // CbcHeuristicRinsFixtureDump_H
