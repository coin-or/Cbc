/**
 * Fixture dumping for CglGomory experiments.
 *
 * @file CbcGomoryFixtureDump.hpp
 * @brief dump the exact inputs of a CglGomory call, for offline replay
 *
 * Same motivation as CbcClqFixtureDump.hpp, CbcZeroHalfFixtureDump.hpp and
 * CbcProbingFixtureDump.hpp: the call worth optimizing only happens deep inside a
 * solve, after preprocessing and the root LP, so experimenting on it costs a full
 * CBC run per measurement. This captures everything the call depends on:
 *
 *   <name>.gomory.mps.gz   the problem as it stands at the call
 *   <name>.gomory.bas      the LP basis -- for Gomory this IS the input, see below
 *   <name>.gomory.sol      the LP solution, for cross-checking the warm start
 *   <name>.gomory.ctype    which columns are integer (MPS cannot say: see below)
 *   <name>.gomory.meta     key/value provenance, read back by the benchmark
 *
 * THE BASIS IS NOT PROVENANCE HERE, IT IS THE ALGORITHM'S INPUT. For the other
 * three generators the `.bas` is a convenience that gets the replay onto the
 * captured vertex cheaply; a cold solve landing on a different optimal vertex
 * would change the experiment but the generator would still be doing the same
 * kind of work. Gomory is different in kind: `CglGomory::needsOptimalBasis()`
 * returns true, and generateCuts' first act is to build a `CoinFactorization` from
 * the basis and BTRAN through it -- the cut *is* a row of the simplex tableau at
 * that basis. Two different optimal bases of the same LP (a degenerate vertex has
 * many) give entirely different tableau rows and therefore different cuts. So a
 * fixture whose LP is not proven optimal is not a weaker Gomory fixture, it is not
 * a Gomory fixture at all, and this dump refuses to write one -- unlike the
 * Probing dump, which writes `lpOptimal 0` and lets the replay solve cold.
 *
 * WHY THE PAYLOAD IS AGAIN JUST THE GENERIC FIVE. Checked by reading
 * CglGomory::refreshSolver() (CglGomory.cpp:2084-2104) rather than assumed: it
 * does nothing but re-take the solver pointer. Every array the workhorse touches
 * -- the column and row copies of the matrix, colsol, the bounds, intVar[] and the
 * warm-start basis -- is passed in as an argument on each call and rebuilt from
 * the solver by the Osi-level entry point. There is nothing cached, so there is no
 * analogue of the clique dump's conflict-graph staleness trap and no
 * `--rebuild-<thing>` flag to expose. Reconstructing from the `.mps` plus the
 * `.bas` is exactly what CBC itself does.
 *
 * The one piece of persistent generator state is `numberTimesStalled_`, and it is
 * dead on this path: the only code that reads it is commented out at
 * CglGomory.cpp:855 and :864. `dynamicLimitInTree_` is set from limit_ but only
 * consulted when `info.inTree`, which is false at the root.
 *
 * WHY NOT REUSE THE probingFixtures OR zerohalfFixtures POPULATION. The five
 * *files* are generator-agnostic and interchangeable; the *set of instances* is
 * not, and Gomory's precondition is the narrowest of the four. It needs a column
 * that is simultaneously **basic**, **integer**, and **fractional by more than
 * `away`** -- see cbcGomoryFixtureCount. Each of the three is a real filter:
 *
 *   - An instance whose root LP is already integral has zero Gomory candidates
 *     and the main loop body never executes, yet it is a perfectly good Probing
 *     fixture (Probing gates on movable bounds, not on fractionality) and a
 *     perfectly good ZeroHalf one.
 *   - A fixed integer column is never basic, so a model of mostly-fixed integers
 *     is a busy Probing case and an empty Gomory one.
 *   - Conversely a general-integer model with fractional coefficients has no row
 *     ZeroHalf will look at, and is a fine Gomory case.
 *
 * So the gate here is Gomory's own, and the population will overlap the other
 * three without matching any of them.
 *
 * The `.ctype` sidecar exists because MPS cannot express "fixed and integer": a
 * column with lb == ub takes writeMps's " FX " branch, which has no integer form,
 * and reads back continuous. Gomory reads `intVar[]` in three separate places, so
 * a lost marker changes the experiment three times over:
 *
 *   - candidate selection (CglGomory.cpp:892) -- a lost marker drops a candidate;
 *   - the coefficient derivation (:1078-1104) -- an integer column contributes
 *     `above_integer(value)` with the mixed-integer Gomory tightening, a
 *     continuous one contributes `value` or `-ratio*value` and also bumps
 *     `numberNonInteger`, which feeds the rhs relaxation ladder at :1606-1638. So
 *     a lost marker does not merely weaken one term, it changes the cut's rhs;
 *   - the integer-slack test (:780) -- `rowType |= 4` requires every entry in the
 *     row to be an integer column, so one lost marker can silence a whole row's
 *     slack contribution.
 *
 * A fixed integer column is admittedly not a *candidate* (it cannot be basic), but
 * it is very much a *coefficient*, which is why this matters for Gomory even
 * though the candidate count would not notice.
 *
 * WHAT THE .meta RECORDS ABOUT THE CALL, AND WHY EACH FIELD IS LOAD-BEARING. A
 * bench that calls generateCuts with a default-constructed CglGomory and a
 * default CglTreeInfo measures a configuration CBC never runs, and in Gomory's
 * case the gap is large. Every field below was derived by reading, and the
 * derivations are recorded because several are non-obvious:
 *
 *   - `cbcAwayAtRoot 0.005` against the constructor's 0.05. The workhorse uses
 *     `away = info.inTree ? away_ : min(away_, awayAtRoot_)`, so at the root the
 *     effective value is 0.005 -- a tenfold wider candidate window than the
 *     library default. This is the single biggest difference between CBC's
 *     configuration and CglGomory's, and it acts directly on the loop count.
 *   - `cbcLimitAtRoot` is 1000, or 2000 when the model has more than 5000
 *     columns (CbcSolverCutSetup.cpp:161-167). "The model" is the *preprocessed*
 *     one: configureCutGenerators runs from babConfigureSearchModel, which is
 *     after preprocess(), so babModel.getNumCols() is the count this dump sees.
 *     Recorded as a resolved number rather than a rule, so the bench need not
 *     re-derive it.
 *   - `cbcLimit 50` is the in-tree cut-length limit. Irrelevant at the root but
 *     recorded, because a bench sweeping `--in-tree` needs it to stay CBC's.
 *   - `infoOptions 0`. This one took the most reading and matters the most, so the
 *     whole derivation is here. CbcCutGenerator.cpp:309-353 builds `info.options`
 *     bit by bit for a non-Probing generator, and at CBC's first root Gomory call
 *     *every* bit is clear:
 *       bit 8 (globalCutsAtRoot) and bit 16 (globalCuts) -- the setters fire only
 *         when howOften < -900 / < -1900 (CbcCutGenerator.cpp:60-100), and
 *         Gomory's howOften is translate[CGIfMove] == -98;
 *       bit 32 (ineffectualCuts) -- set only from strongCuts (CbcModel.cpp:9543),
 *         which needs a populated cut_obj history; it is all -COIN_DBL_MAX on the
 *         first pass, so the test fails;
 *       bit 128 -- set when fullScan < 0, but at the root numberNodes_ == 0 gives
 *         fullScan == 2 (CbcModel.cpp:10956-10962);
 *       bit 256 -- needs moreSpecialOptions() & 16384, which nothing in Cbc/src
 *         sets;
 *       bit 512 -- needs a parentModel;
 *       bit 1024 -- needs must-call-again mode.
 *     Also `switches[Gomory] = lagrangeanFlag`, which is 0 by default, and the
 *     apply loop masks it with `& ~16383` before or-ing, so nothing arrives there
 *     either. Four consequences inside the workhorse, each worth knowing before
 *     optimizing anything:
 *       `testFixed = 1.0e-8` rather than -1.0, so fixed columns are excluded from
 *         the dense tableau loop;
 *       `doSorted` (bit 256) is false, which means nTotalEls is COIN_INT_MAX (no
 *         element budget), candidates are NOT sorted by fractionality, and the
 *         MORE_GOMORY_CUTS==1 `secondaryCuts` flush at :1893-1922 never runs --
 *         so the `accurate2` cuts are computed and then silently discarded;
 *       no cut is marked globally valid (:337-347 needs bit 4, 8 or 16);
 *       the unconditional printf at :81 is gated on bit 16 and stays quiet.
 *   - `infoPass 0`. The tolerance regime at :844-853 keys off it: pass 0 sets
 *     tolerance1 = 1.0, which makes the accuracy test at :1238-1259 *absolute*
 *     (1e-7) rather than relative, and tolerance2 = 1e-2, so acceptance needs
 *     `sum > rhs + 1e-2*away`. Later passes use a relative test. This is a
 *     different acceptance rule, not a tighter one.
 *   - `cbcGomoryType 0` and no original solver. Together these make the whole
 *     CGL_HAS_CLP_GOMORY prologue at :100-275 dead code on CBC's default root
 *     path: `clpSolver` comes out NULL and `useSolver` stays `&si`. A bench that
 *     called passInOriginalSolver would measure the lagrangean variants CBC only
 *     uses under -lagomory.
 *   - `maxSeconds` is the dumping run's budget, provenance only. CbcCutGenerator
 *     does call setMaxSeconds(remaining) on this generator (:349), but grepping
 *     CglGomory.cpp for maxSeconds/CoinCpuTime/getTimeOfDay finds nothing -- the
 *     setter lands on the base class and the value is never read. So the root
 *     Gomory call has no internal deadline and a replay is reproducible.
 *
 * Entirely behind CBC_DUMP_GOMORY_FIXTURE and off by default. Header-only static
 * functions, so no Makefile.am/Makefile.in changes are needed.
 *
 * Environment:
 *   CBC_GOMORY_FIXTURE_DIR    output directory
 *                             (default ~/instances/miplib/2017+spp/gomoryFixtures)
 *   CBC_GOMORY_FIXTURE_NAME   instance base name, when the MPS problem name is
 *                             unhelpful (drivers normally set this)
 **/

#ifndef CbcGomoryFixtureDump_H
#define CbcGomoryFixtureDump_H

#ifdef CBC_DUMP_GOMORY_FIXTURE

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
static void cbcGomoryFixtureMkdirP(const std::string &path)
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
static std::string cbcGomoryFixtureDir()
{
  const char *env = getenv("CBC_GOMORY_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/miplib/2017+spp/gomoryFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback, but
/// drivers should set CBC_GOMORY_FIXTURE_NAME: several instances share a problem
/// name, and some carry none at all.
static std::string cbcGomoryFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_GOMORY_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

/**
 * Write the MPS. formatType 2 is IEEE hex, which round-trips every double exactly;
 * plain MPS loses bits. For Gomory that is not a fussy detail -- the cut
 * coefficients come out of a factorization of this matrix, and the acceptance test
 * at :1238-1259 compares two independently accumulated sums against an absolute
 * 1e-7. A matrix perturbed in its last bits can flip a cut from accepted to
 * rejected. numberAcross 1 for the same reason.
 *
 * Empty columns are the one wrinkle: writeMps drops a column with no matrix entry,
 * which shifts every later column index and so invalidates the .ctype sidecar, the
 * .bas and the .sol. A single redundant final row covering those columns keeps them
 * at their original index; `paddedColumns` in the .meta says how many there were,
 * and the loader deletes the extra row.
 */
static int cbcGomoryFixtureWriteMps(const OsiSolverInterface *si,
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
    // activity can take, so it is redundant by construction. Deliberately not
    // truly free -- a doubly-infinite row would be dropped by some generators and
    // would take the column indices with it. Note the row is basic in the .bas
    // written next (it is slack in the captured basis by construction, since it
    // did not exist then and the loader deletes it before reading the basis), so
    // it cannot perturb the factorization either.
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
  // success*, leaving the replay to start cold -- which for Gomory means no cuts
  // at all rather than merely a different vertex, so this is load-bearing twice
  // over.
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
static bool cbcGomoryFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
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

/// Write the LP solution in the same shape as preProcessedInstances/*.sol: an
/// "=obj=" line, then one line per nonzero as "<index> <name> <value>".
static bool cbcGomoryFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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
 * `isContinuous()` is the source of truth rather than `getColType()`, which derives
 * Binary/GeneralInteger from the current bounds.
 */
static bool cbcGomoryFixtureWriteColTypes(const OsiSolverInterface *si,
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

/// CglGomory's own fractionality test (CglGomory.cpp:363-369), reproduced so the
/// candidate count is the generator's and not an approximation of it. The 1e-9
/// relative snap matters: a column at 3.9999999997 is *not* a candidate.
static double cbcGomoryFixtureAboveInteger(double value)
{
  const double value2 = floor(value);
  const double value3 = floor(value + 0.5);
  if (fabs(value3 - value) < 1.0e-9 * (fabs(value3) + 1.0))
    return 0.0;
  return value - value2;
}

/**
 * Count what CglGomory would actually work on, mirroring its own tests.
 *
 * `gomoryCandidates` is the acceptance test that decides whether the generator can
 * do anything at all, and it is the skip gate. CglGomory.cpp:889-901 keeps column
 * iColumn when it is **basic**, **integer**, and `above_integer(colsol)` lies
 * strictly inside (away, 1-away). Zero candidates means the main loop body never
 * executes: no BTRAN, no tableau row, no cut. Nothing to measure. Note this is a
 * strictly narrower gate than either Probing's (movable integer bounds) or
 * ZeroHalf's (a row of integral coefficients on integer columns), which is why
 * this population cannot be inherited from theirs.
 *
 * `basicIntegers` is reported alongside so a skip can be read: `basicIntegers 300,
 * gomoryCandidates 0` means the root LP is integral on every basic integer column
 * -- a genuinely cut-free instance -- whereas `basicIntegers 0` means the integer
 * columns are all nonbasic or fixed, a different situation with the same verdict.
 *
 * `nonbasicMovable` is the width of the dense loop at CglGomory.cpp:1037-1115,
 * which reproduces its own guard exactly: `columnIsBasic[j] < 0 && colUpper[j] >
 * colLower[j] + testFixed` with testFixed = 1e-8, the value that applies when
 * options bit 16 is clear -- i.e. at CBC's root. That loop runs once per candidate
 * and dots the BTRAN result against every element of every such column, so
 * `gomoryCandidates * nonbasicMovableNz` is the leading term in the separation
 * cost and is the number to rank the slow tail by. It is written to the .meta as
 * `denseLoopWork` so a scan script need not multiply.
 *
 * `activeSlackRows` / `integerSlackRows` reproduce rowType (CglGomory.cpp:747-795),
 * which decides the slack part of the cut at :1131-1176. A row contributes only
 * when it is nonbasic and not an equality; the `|4` subset (integral rhs, all
 * entries integer columns with integral coefficients) additionally gets the
 * integer-slack tightening. Both need the row activity at colsol, so it is computed
 * here the same way the generator does at :728-738.
 */
static void cbcGomoryFixtureCount(const OsiSolverInterface *si,
  const CoinWarmStartBasis *ws, double away, int &gomoryCandidates,
  int &basicIntegers, int &nonbasicMovable, double &nonbasicMovableNz,
  int &activeSlackRows, int &integerSlackRows)
{
  gomoryCandidates = basicIntegers = nonbasicMovable = 0;
  activeSlackRows = integerSlackRows = 0;
  nonbasicMovableNz = 0.0;

  if (!ws)
    return;

  const int numberColumns = si->getNumCols();
  const int numberRows = si->getNumRows();
  const double *colLower = si->getColLower();
  const double *colUpper = si->getColUpper();
  const double *rowLower = si->getRowLower();
  const double *rowUpper = si->getRowUpper();
  const double *colsol = si->getColSolution();
  if (!colsol)
    return;

  // A basis written for a different model would index out of range; refuse rather
  // than count nonsense.
  if (ws->getNumStructural() != numberColumns || ws->getNumArtificial() != numberRows)
    return;

  // testFixed is 1.0e-8 when options bit 16 is clear, which is CBC's root case.
  const double testFixed = 1.0e-8;

  const CoinPackedMatrix *byCol = si->getMatrixByCol();
  const int *columnLength = byCol ? byCol->getVectorLengths() : NULL;

  for (int j = 0; j < numberColumns; ++j) {
    const bool basic
      = ws->getStructStatus(j) == CoinWarmStartBasis::basic;
    const bool isInt = !si->isContinuous(j);
    if (basic) {
      if (isInt) {
        ++basicIntegers;
        const double frac = cbcGomoryFixtureAboveInteger(colsol[j]);
        if (frac < 1.0 - away && frac > away)
          ++gomoryCandidates;
      }
    } else if (colUpper[j] > colLower[j] + testFixed) {
      ++nonbasicMovable;
      if (columnLength)
        nonbasicMovableNz += (double)columnLength[j];
    }
  }

  if (!byCol)
    return;

  // Row activity at colsol, exactly as CglGomory.cpp:728-738 builds it.
  const int *row = byCol->getIndices();
  const CoinBigIndex *columnStart = byCol->getVectorStarts();
  const double *columnElements = byCol->getElements();
  std::vector< double > rowActivity(numberRows, 0.0);
  for (int j = 0; j < numberColumns; ++j) {
    const double value = colsol[j];
    for (CoinBigIndex k = columnStart[j]; k < columnStart[j] + columnLength[j]; ++k)
      rowActivity[row[k]] += columnElements[k] * value;
  }

  const CoinPackedMatrix *byRow = si->getMatrixByRow();
  const int *column = byRow ? byRow->getIndices() : NULL;
  const CoinBigIndex *rowStart = byRow ? byRow->getVectorStarts() : NULL;
  const int *rowLength = byRow ? byRow->getVectorLengths() : NULL;
  const double *rowElements = byRow ? byRow->getElements() : NULL;

  for (int i = 0; i < numberRows; ++i) {
    const bool rowBasic = ws->getArtifStatus(i) == CoinWarmStartBasis::basic;
    if (rowBasic || rowUpper[i] <= rowLower[i] + 1.0e-7)
      continue; // equality or basic: rowType 0, contributes nothing
    double rhs = 0.0;
    if (rowActivity[i] >= rowUpper[i] - 1.0e-7) {
      rhs = rowUpper[i];
    } else if (rowActivity[i] <= rowLower[i] + 1.0e-7) {
      rhs = rowLower[i];
    } else {
      // The generator's "probably large rhs" branch: it sets rowType but then
      // `continue`s past the integer-slack test, so the row is active without the
      // |4 bit.
      ++activeSlackRows;
      continue;
    }
    ++activeSlackRows;
    if (!byRow || cbcGomoryFixtureAboveInteger(rhs) >= 1.0e-10)
      continue;
    bool allInteger = true;
    for (CoinBigIndex k = rowStart[i]; k < rowStart[i] + rowLength[i]; ++k) {
      const int j = column[k];
      if (si->isContinuous(j)
        || cbcGomoryFixtureAboveInteger(rowElements[k]) > 1.0e-10) {
        allInteger = false;
        break;
      }
    }
    if (allInteger)
      ++integerSlackRows;
  }
}

/**
 * Write the provenance file.
 *
 * Beyond the usual model shape this records the *call shape* and CBC's knob
 * values; see the header comment for the derivation of each, and for why a bench
 * that guesses them measures a call CBC never makes.
 */
static bool cbcGomoryFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  int paddedColumns, int integerColumns, int gomoryCandidates, int basicIntegers,
  int nonbasicMovable, double nonbasicMovableNz, int activeSlackRows,
  int integerSlackRows, int pass, int options, int inTree, int formulationRows,
  int limitAtRoot, double awayAtRoot, double cutoff, double maxSeconds,
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
  // The dumping run's budget, NOT a limit the Gomory call honoured: CbcCutGenerator
  // does call setMaxSeconds on this generator, but CglGomory never reads it (no
  // maxSeconds/CoinCpuTime/getTimeOfDay reference anywhere in CglGomory.cpp), so
  // the call has no internal deadline and a replay is reproducible.
  fprintf(fp, "maxSeconds %.6f\n", maxSeconds);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index. `rows` above is
  // the captured count, so a consumer loading rows+1 deletes the last row.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());
  // What CglGomory would look at: see cbcGomoryFixtureCount.
  fprintf(fp, "gomoryCandidates %d\n", gomoryCandidates);
  fprintf(fp, "basicIntegers %d\n", basicIntegers);
  fprintf(fp, "nonbasicMovable %d\n", nonbasicMovable);
  fprintf(fp, "nonbasicMovableNz %.0f\n", nonbasicMovableNz);
  fprintf(fp, "activeSlackRows %d\n", activeSlackRows);
  fprintf(fp, "integerSlackRows %d\n", integerSlackRows);
  // The leading term in separation cost: the dense loop at CglGomory.cpp:1037-1115
  // runs once per candidate over every element of every nonbasic movable column.
  // Precomputed so a ranking script does not have to multiply, and as a double
  // because the product overflows int on the larger instances.
  fprintf(fp, "denseLoopWork %.0f\n", (double)gomoryCandidates * nonbasicMovableNz);
  // The call shape CBC used, so a bench can reproduce it rather than guess.
  fprintf(fp, "infoPass %d\n", pass);
  fprintf(fp, "infoOptions %d\n", options);
  fprintf(fp, "infoInTree %d\n", inTree);
  fprintf(fp, "infoFormulationRows %d\n", formulationRows);
  fprintf(fp, "cutoff %.15g\n", cutoff);
  // CBC's Gomory configuration (CbcSolverCutSetup.cpp:132-188), which differs from
  // CglGomory's constructor defaults most sharply on awayAtRoot: 0.005 against the
  // library's 0.05, a tenfold wider candidate window. limitAtRoot is resolved here
  // rather than left as a rule, because the >5000 test is on the *preprocessed*
  // column count, which only this call site knows.
  fprintf(fp, "cbcLimitAtRoot %d\n", limitAtRoot);
  fprintf(fp, "cbcLimit 50\n");
  fprintf(fp, "cbcAwayAtRoot %.15g\n", awayAtRoot);
  fprintf(fp, "cbcAway 0.05\n");
  fprintf(fp, "cbcGomoryType 0\n");
  fprintf(fp, "cbcAlternateFactorization 0\n");
  fprintf(fp, "cbcConditionNumberMultiplier 1e-18\n");
  fprintf(fp, "cbcLargestFactorMultiplier 1e-13\n");
  fclose(fp);
  return true;
}

/**
 * Dump one fixture. Returns true if files were written.
 *
 * Two skip criteria, and both are Gomory's own rather than inherited:
 *
 *   - no proven-optimal LP. Unlike the Probing dump, which writes `lpOptimal 0`
 *     and lets the replay solve cold, a Gomory fixture without a basis is useless:
 *     needsOptimalBasis() is true and the cut is a row of the tableau at that
 *     basis. A cold re-solve would land on some optimal basis, but not the one CBC
 *     used, and on a degenerate vertex that is a different set of cuts.
 *   - no candidate: no column is simultaneously basic, integer, and fractional by
 *     more than away = min(away_, awayAtRoot_) = 0.005. The main loop body then
 *     never executes and there is nothing to measure.
 */
static bool cbcDumpGomoryFixture(OsiSolverInterface *si, const char *tag = "gomory",
  int pass = 0, int options = 0, int inTree = 0, int formulationRows = -1,
  double maxSeconds = 0.0)
{
  if (!si)
    return false;

  const std::string name = cbcGomoryFixtureName(si);

  if (si->getNumCols() == 0 || si->getNumRows() == 0) {
    printf("[gomoryfixture] %s.%s: SKIP (empty model: %d rows, %d cols)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols());
    fflush(stdout);
    return false;
  }

  if (!si->isProvenOptimal()) {
    printf("[gomoryfixture] %s.%s: SKIP (LP not proven optimal; CglGomory's "
           "needsOptimalBasis() is true and the cut is a tableau row at that "
           "basis, so a fixture without one is not replayable)\n",
      name.c_str(), tag);
    fflush(stdout);
    return false;
  }

  const CoinWarmStartBasis *ws
    = dynamic_cast< const CoinWarmStartBasis * >(si->getWarmStart());
  if (!ws) {
    printf("[gomoryfixture] %s.%s: SKIP (no CoinWarmStartBasis from the solver)\n",
      name.c_str(), tag);
    fflush(stdout);
    return false;
  }

  // CBC's effective root value: the workhorse uses min(away_, awayAtRoot_) when
  // !info.inTree, and CBC sets away_ = 0.05 (constructor) with awayAtRoot_ = 0.005.
  const double awayAtRoot = 0.005;
  const double away = awayAtRoot < 0.05 ? awayAtRoot : 0.05;

  int gomoryCandidates = 0, basicIntegers = 0, nonbasicMovable = 0;
  int activeSlackRows = 0, integerSlackRows = 0;
  double nonbasicMovableNz = 0.0;
  cbcGomoryFixtureCount(si, ws, away, gomoryCandidates, basicIntegers,
    nonbasicMovable, nonbasicMovableNz, activeSlackRows, integerSlackRows);
  delete ws;

  if (gomoryCandidates == 0) {
    printf("[gomoryfixture] %s.%s: SKIP (no Gomory candidate; %d rows, %d cols, "
           "%d basic integer column(s) -- none has above_integer(colsol) inside "
           "(%g, %g), so the main loop is empty)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols(), basicIntegers,
      away, 1.0 - away);
    fflush(stdout);
    return false;
  }

  double cutoff = COIN_DBL_MAX;
  si->getDblParam(OsiDualObjectiveLimit, cutoff);

  // CbcSolverCutSetup.cpp:161-167, on the *preprocessed* column count -- which is
  // what this solver holds, since configureCutGenerators runs after preprocess().
  const int limitAtRoot = si->getNumCols() > 5000 ? 2000 : 1000;

  const std::string dir = cbcGomoryFixtureDir();
  cbcGomoryFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcGomoryFixtureWriteMps(si, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk
    = cbcGomoryFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const bool metaOk = cbcGomoryFixtureWriteMeta(si, tag, paddedColumns,
    integerColumns, gomoryCandidates, basicIntegers, nonbasicMovable,
    nonbasicMovableNz, activeSlackRows, integerSlackRows, pass, options, inTree,
    formulationRows >= 0 ? formulationRows : si->getNumRows(), limitAtRoot,
    awayAtRoot, cutoff, maxSeconds, stem + ".meta");
  const bool basOk = cbcGomoryFixtureWriteBasis(si, stem + ".bas");
  const bool solOk = cbcGomoryFixtureWriteSol(si, stem + ".sol");

  printf("[gomoryfixture] %s.%s: %s rows=%d cols=%d int=%d padded=%d cand=%d "
         "basicInt=%d nbMovable=%d nbMovableNz=%.0f activeSlack=%d intSlack=%d "
         "denseWork=%.0f limitAtRoot=%d obj=%.15g pass=%d options=%d "
         "(mps=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, paddedColumns,
    gomoryCandidates, basicIntegers, nonbasicMovable, nonbasicMovableNz,
    activeSlackRows, integerSlackRows,
    (double)gomoryCandidates * nonbasicMovableNz, limitAtRoot,
    si->getObjValue(), pass, options,
    mpsRc, (int)ctOk, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  // The basis is required, not optional, so unlike the Probing dump it counts
  // toward success.
  return mpsRc == 0 && ctOk && metaOk && basOk;
}

#endif // CBC_DUMP_GOMORY_FIXTURE
#endif // CbcGomoryFixtureDump_H
