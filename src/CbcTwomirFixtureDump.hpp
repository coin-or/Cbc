/**
 * Fixture dumping for CglTwomir experiments.
 *
 * @file CbcTwomirFixtureDump.hpp
 * @brief dump the exact inputs of a CglTwomir call, for offline replay
 *
 * Same motivation as CbcClqFixtureDump.hpp, CbcZeroHalfFixtureDump.hpp,
 * CbcProbingFixtureDump.hpp and CbcGomoryFixtureDump.hpp: the call worth
 * optimizing only happens deep inside a solve, after preprocessing and the root
 * LP, so experimenting on it costs a full CBC run per measurement. This captures
 * everything the call depends on:
 *
 *   <name>.twomir.mps.gz   the problem as it stands at the call
 *   <name>.twomir.bas      the LP basis -- for Twomir this IS the input, see below
 *   <name>.twomir.sol      the LP solution, for cross-checking the warm start
 *   <name>.twomir.ctype    which columns are integer (MPS cannot say: see below)
 *   <name>.twomir.meta     key/value provenance, read back by the benchmark
 *
 * THE BASIS IS THE ALGORITHM'S INPUT, NOT PROVENANCE -- as for Gomory, and by a
 * different route. `CglTwomir::needsOptimalBasis()` returns true
 * (CglTwomir.cpp:2241-2246), and `DGG_generateTabRowCuts` builds its *own*
 * `CoinFactorization` from that basis (:1569) and BTRANs a unit vector through it
 * (:1044) once per candidate; the tableau row is then formed by dotting the
 * resulting B^-1 row against the column-major matrix (:1056-1067). So Twomir does
 * not go through `getBInvARow`/`getBInvRow` at all, but it is just as much a
 * function of the specific basis: two different optimal bases of the same LP (a
 * degenerate vertex has many) give different tableau rows and therefore different
 * cuts. A fixture whose LP is not proven optimal is not a weaker Twomir fixture,
 * it is not a Twomir fixture at all, and this dump refuses to write one.
 *
 * WHY THE PAYLOAD IS AGAIN JUST THE GENERIC FIVE. Checked by reading
 * `CglTwomir::refreshSolver()` (CglTwomir.cpp:2269-2279) rather than assumed. It
 * is *not* the no-op CglGomory's is -- it clones the solver into
 * `originalSolver_` -- but that path is entirely conditional on
 * `originalSolver_` already being non-NULL, which only `passInOriginalSolver()`
 * makes true, and CBC only calls that under `-latwomir`
 * (CbcSolverCutSetup.cpp:402). On the default root path `originalSolver_` is NULL,
 * `clpSolver` comes out NULL at :99, and the whole CGL_HAS_CLP_TWOMIR prologue
 * (:97-281) is dead: `useSolver` stays `&si`. Everything the workhorse touches
 * -- `DGG_getData`'s bounds/solution/reduced-cost/dual arrays, both matrix
 * copies, `getColType()` -- is re-read from the solver on every call, and
 * `DGG_freeData` releases all of it. There is nothing cached across calls, so
 * there is no analogue of the clique dump's conflict-graph staleness trap.
 *
 * The one piece of persistent generator state is `randomNumberGenerator_`
 * (CglTwomir.cpp:320), passed to the formulation stage by non-const reference; its
 * `seed_` is `mutable` and each `randomDouble()` advances an LCG. It is
 * constructed with seed 987654321 at :569 -- but note that under
 * `COIN_OWN_RANDOM_32`, which is defined unconditionally at
 * CoinHelperFunctions.hpp:892, `CoinThreadRandom(unsigned)` honours the seed, so
 * 987654321 it is, and `info.randomNumberGenerator` (which CbcCutGenerator does
 * fill in at :323) is *ignored* by CglTwomir. CBC constructs the generator once
 * and calls it every pass, so pass 1 sees the pristine stream; a bench that
 * replays a fixture with a fresh CglTwomir therefore reproduces CBC's pass-1
 * draws exactly, and a multi-round bench does not reproduce pass 2's.
 *
 * WHY THE GOMORY OR PROBING POPULATION CANNOT BE REUSED. The five *files* are
 * generator-agnostic and interchangeable; the *set of instances* is not. Twomir's
 * preconditions are neither a subset nor a superset of Gomory's:
 *
 *   - A FREE COLUMN DISQUALIFIES THE WHOLE INSTANCE. CglTwomir.cpp:109-112 counts
 *     columns with `colLower < -1e20 && colUpper > 1e20`, and if there is even one
 *     it `return`s at :118 having produced nothing at all. No other generator here
 *     does this. `freeColumns` is recorded and is its own skip arm, because it is
 *     the single biggest reason this population differs from Gomory's.
 *   - CONVERSELY, `twomirCandidates == 0` IS NOT A VALID SKIP, unlike Gomory's
 *     `gomoryCandidates == 0`. Twomir has two independent cut sources: the tableau
 *     stage (`do_tab_`, :314) needs a basic fractional integer column, but the
 *     formulation stage (`do_form_`, :317) walks `info.formulation_rows` rows of
 *     the *original formulation* and needs no fractional basic integer whatever.
 *     A Gomory-shaped skip here would silently drop the entire formulation-cut
 *     population, so the gate is `twomirCandidates == 0 && formIntRows == 0`.
 *   - Twomir does NOT exclude fixed nonbasic columns from its dense loop.
 *     CglGomory's equivalent loop guards on `colUpper[j] > colLower[j] + 1e-8`;
 *     CglTwomir.cpp:1057 guards only on `!DGG_isBasic(data, j) || j == index`. So
 *     an instance of mostly-fixed columns is a *cheap* Gomory case and an
 *     expensive Twomir one, and `nonbasicStructuralNz` below is deliberately not
 *     filtered by fixedness.
 *
 * The `.ctype` sidecar exists because MPS cannot express "fixed and integer": a
 * column with lb == ub takes writeMps's " FX " branch, which has no integer form,
 * and reads back continuous. For Twomir a lost marker changes the experiment in
 * four places, one of which is worse than anything in the Gomory case:
 *
 *   - `DGG_getData` (:802-808) sets the integer flag *and tightens the stored
 *     bounds* with `ceil(colLower)` / `floor(colUpper)`. A column read back as
 *     continuous keeps untightened bounds, and `data->lb`/`data->ub` feed
 *     `DGG_transformConstraint`'s complement decision (:1345-1358), so the cut is
 *     built around a different origin. This is not a weakening, it is a different
 *     cut.
 *   - candidate selection (:1575) -- a lost marker drops a candidate;
 *   - the row-integrality test (:918) -- `DGG_setIsInteger` for a *slack* requires
 *     every entry in the row to be an integer column, so one lost marker can
 *     silence a whole row's integer-slack status, which in turn removes the row
 *     from `formIntRows` and changes `isint[]` in the formulation stage;
 *   - `intVar[iColumn]` in the singleton branch (:429) decides whether a
 *     one-element cut becomes an integer-rounded `OsiColCut` or is discarded.
 *
 * WHAT THE .meta RECORDS ABOUT THE CALL, AND WHY EACH FIELD IS LOAD-BEARING. A
 * bench that calls generateCuts with a default-constructed CglTwomir and a default
 * CglTreeInfo measures a configuration CBC never runs. For Twomir the gap runs in
 * the *opposite* direction from Gomory's: the library default `awayAtRoot_` is
 * 0.0005 (:570) and CBC widens it to 0.005, a tenfold *narrower* candidate window.
 * Derivations, several non-obvious:
 *
 *   - `cbcAwayAtRoot 0.005` / `cbcAway 0.01` come from
 *     CbcSolverCutSetup.cpp:365-367, the `#ifdef MORE_CUTS` branch. That branch is
 *     live for a thoroughly fragile reason: `#define MORE_CUTS` sits at :139,
 *     inside the *Gomory* configuration block, and is never `#undef`'d. Deleting
 *     it there would silently change Twomir's root window to 0.01.
 *   - `cbcMaxElements 250` (:362) and `cbcMaxElementsRoot 50000` -- the latter is
 *     the constructor's value, because **`setMaxElementsRoot` is never called
 *     anywhere in Cbc**. That matters because :296 selects between them by
 *     `info.inTree`, so at the root the 250 is *not* the cap.
 *   - `effMaxElements` is what `max_elements` actually becomes, and at pass 0 it is
 *     neither 250 nor 50000: :301-302 overwrites it with
 *     `useSolver->getNumCols()` when `!info.pass || (info.options & 32)`. It gates
 *     the emit loop at :335 (`cut->nz < max_elements`), so it decides how many of
 *     the generated cuts survive to `cs`.
 *   - `infoOptions 0`, re-derived for Twomir rather than inherited from the Gomory
 *     dump. Every bit CbcCutGenerator.cpp:309-353 can set is clear at the first
 *     root pass: Twomir's howOften is `translate[CGIfMove] == -98`, which is not
 *     `< -900` or `< -1900`, so bits 8 and 16 stay clear; `strongCuts` is false on
 *     pass 1 because the cut_obj history is still -COIN_DBL_MAX (bit 32); fullScan
 *     is 2 at the root (bit 128); nothing in Cbc/src sets
 *     `moreSpecialOptions() & 16384` (bit 256); no parentModel (bit 512); not
 *     must-call-again (bit 1024). Bit 32 matters *more* here than for Gomory: it
 *     is the second half of the `max_elements = getNumCols()` test at :301, so on a
 *     later pass with strongCuts set the cap jumps back up rather than down.
 *     `switches[Twomir] = 1 | lagrangeanFlag` (:383) and lagrangeanFlag is 0 by
 *     default; the apply loop masks with `& ~16383` before or-ing, so the 1 (which
 *     is `setSwitchOffIfLessThan(1)`, see below) does not reach info.options.
 *   - `infoPass 0` and `infoLevel 0`. Both are load-bearing and in different ways.
 *     `info.pass` selects the `max_elements` override just described. `info.level`
 *     and `info.pass` *together* gate the tableau stage: :314 runs it only when
 *     `do_tab_ && info.level < 1 && info.pass < 6`. So pass 6 onward is a
 *     structurally different call -- the formulation stage alone -- not a cheaper
 *     version of this one.
 *   - `infoFormulationRows` is `model_->numberRowsAtContinuous()`
 *     (CbcCutGenerator.cpp:312), and it caps the formulation stage's outer loop at
 *     :1617 (`num_rows = min(data->nrow, nrows)`). It is a call-shape input the
 *     generator cannot derive, and note that `setFormulationRows()` does *not*
 *     supply it: `form_nrows_` is written by that setter and never read.
 *   - `cbcSwitchOffIfLessThan 1`. `switches[...] = 1` at :383 becomes
 *     `setSwitchOffIfLessThan(1)`, and CbcCutGenerator.cpp:1584-1592 then disables
 *     TwoMirCuts *for the whole solve* if this call returns fewer than 1 row cut.
 *     Two consequences: the population this dump defines is exactly the population
 *     the kill switch judges, and a scan of CBC's generator table must distinguish
 *     "TwoMirCuts row absent" (never added) from "row present, Next run =
 *     disabled" (added, fired once, killed).
 *   - `maxSeconds` is the dumping run's budget, provenance only. CbcCutGenerator
 *     does call `setMaxSeconds(remaining)` on this generator, but grepping
 *     CglTwomir.{cpp,hpp} for maxSeconds / CoinCpuTime / CoinWallclockTime /
 *     getTimeOfDay / clock finds *nothing* -- the setter lands on the base class
 *     and the value is never read. So the root Twomir call has no internal deadline
 *     and a replay is reproducible.
 *
 * Entirely behind CBC_DUMP_TWOMIR_FIXTURE and off by default. Header-only static
 * functions, so no Makefile.am/Makefile.in changes are needed.
 *
 * Environment:
 *   CBC_TWOMIR_FIXTURE_DIR    output directory
 *                             (default ~/instances/miplib/2017+spp/twomirFixtures)
 *   CBC_TWOMIR_FIXTURE_NAME   instance base name, when the MPS problem name is
 *                             unhelpful (drivers normally set this)
 **/

#ifndef CbcTwomirFixtureDump_H
#define CbcTwomirFixtureDump_H

#ifdef CBC_DUMP_TWOMIR_FIXTURE

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
static void cbcTwomirFixtureMkdirP(const std::string &path)
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
static std::string cbcTwomirFixtureDir()
{
  const char *env = getenv("CBC_TWOMIR_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/miplib/2017+spp/twomirFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback, but
/// drivers should set CBC_TWOMIR_FIXTURE_NAME: several instances share a problem
/// name, and some carry none at all.
static std::string cbcTwomirFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_TWOMIR_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

/**
 * Write the MPS. formatType 2 is IEEE hex, which round-trips every double exactly;
 * plain MPS loses bits. For Twomir that is not a fussy detail -- the tableau row
 * coefficients come out of a factorization of this matrix, and the pipeline that
 * consumes them is threshold-dense: DGG_MIN_TABLEAU_COEFFICIENT (1e-12) decides
 * which coefficients exist at all, DGG_nicefyConstraint rounds near-integers and
 * pads the rhs by up to DGG_NICEFY_MAX_PADDING (1e-6), and DGG_MIN_RHO /
 * DGG_MIN_ALPHA (both 1e-5) decide whether a 2-step MIR is built. A matrix
 * perturbed in its last bits can flip a cut from generated to absent.
 * numberAcross 1 for the same reason.
 *
 * Empty columns are the one wrinkle: writeMps drops a column with no matrix entry,
 * which shifts every later column index and so invalidates the .ctype sidecar, the
 * .bas and the .sol. A single redundant final row covering those columns keeps them
 * at their original index; `paddedColumns` in the .meta says how many there were,
 * and the loader deletes the extra row.
 */
static int cbcTwomirFixtureWriteMps(const OsiSolverInterface *si,
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
    // would take the column indices with it. The row is slack-basic in the .bas
    // written next by construction (it did not exist when the basis was taken, and
    // the loader deletes it before reading the basis), so it cannot perturb the
    // factorization either.
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
  // success*, leaving the replay to start cold -- which for Twomir means a
  // different factorization and so a different cut set, not merely a different
  // vertex, so this is load-bearing twice over.
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
static bool cbcTwomirFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
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
static bool cbcTwomirFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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
 * Binary/GeneralInteger from the current bounds -- and note that CglTwomir reads
 * `getColType()` in three places, one of which (DGG_getData) uses it to *tighten*
 * the bounds it stores, so a stale or lost marker changes the cut rather than
 * merely weakening it. See the header comment.
 */
static bool cbcTwomirFixtureWriteColTypes(const OsiSolverInterface *si,
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

/// CglTwomir's own fractionality function (CglTwomir.hpp:456-459), reproduced so
/// the candidate count is the generator's and not an approximation of it. Note
/// there is deliberately **no** near-integer snap here, unlike CglGomory's
/// above_integer: a value at 3.9999999997 yields 0.9999999997, which the
/// (threshold, 1-threshold) window then rejects at the top end rather than the
/// bottom. The same asymmetry governs the coefficient-integrality tests below,
/// where `frac_part(-3.0) == 0` but `frac_part(-3.0000000001) ~= 1`, so a
/// negative coefficient a hair below an integer is *not* treated as integral.
static double cbcTwomirFixtureFracPart(double value)
{
  return value - floor(value);
}

/// Everything cbcTwomirFixtureCount produces. A struct rather than seventeen out
/// parameters, and the field names are exactly the .meta keys.
struct CbcTwomirFixtureCounts {
  int freeColumns;
  int integerColumns;
  int fixedColumns;
  int basicIntegers;
  int twomirCandidates;
  int nonbasicStructural;
  double nonbasicStructuralNz;
  double tabRowWork;
  double outerLoopWork;
  int basicRows;
  int equalityRows;
  int activeSlackRows;
  int negativeSlackRows;
  int integerSlackRows;
  int formBaseRows;
  int formIntRows;
  int rowNzMax;
};

/**
 * Count what CglTwomir would actually work on, mirroring DGG_getData
 * (CglTwomir.cpp:702-936) and DGG_generateTabRowCuts' candidate test (:1575-1579)
 * rather than approximating them.
 *
 * `freeColumns` IS THE FIRST GATE AND IT IS ALL-OR-NOTHING. CglTwomir.cpp:109-112
 * counts columns with `colLower < -1e20 && colUpper > 1e20`; if the count is
 * nonzero, generateCuts `return`s at :118 with no cuts at all, whatever else the
 * instance looks like. No other generator in this workspace behaves this way, and
 * this is the main reason the Twomir population is not the Gomory population.
 *
 * `twomirCandidates` is the tableau stage's work list: :1575-1579 keeps column k
 * when it is **basic**, **integer**, and `frac_part(x_k)` lies strictly inside
 * (threshold, 1-threshold) with threshold = `data->gomory_threshold` =
 * `awayAtRoot_` = 0.005 at the root (:297, not min()'d with away_, unlike
 * CglGomory). `basicIntegers` is reported alongside so a zero can be read: `
 * basicIntegers 300, twomirCandidates 0` means the root LP is integral on every
 * basic integer column, whereas `basicIntegers 0` means they are all nonbasic.
 *
 * `nonbasicStructural` / `nonbasicStructuralNz` are the width and the element
 * count of the dense loop at :1056-1067, reproducing its guard exactly:
 * `!DGG_isBasic(data, j) || j == index`. **Fixedness is deliberately not
 * filtered**, because Twomir does not filter it -- CglGomory's analogous loop
 * excludes `colUpper[j] <= colLower[j] + 1e-8` and Twomir's does not. That loop
 * runs once per candidate and dots the BTRAN result against every element of every
 * such column, so `tabRowWork = twomirCandidates * nonbasicStructuralNz` is the
 * leading term in tableau-stage cost and is the number to rank the slow tail by.
 *
 * `outerLoopWork = twomirCandidates * cols` is the **floor** that the same loop
 * cannot go below, and it is recorded for a specific reason: any sparsification of
 * the inner element loop must keep the outer `for (j = 0; j < ncol; j++)` sweep
 * dense and in increasing j, because CoinSort_2 at :1099 re-sorts only the slack
 * tail (`tabrow->index + nz1` onward) and the structural head's order is part of
 * the contract its consumers rely on. A fixture where `outerLoopWork` dominates
 * `tabRowWork` therefore *cannot* improve, and without publishing both, such a
 * fixture looks like a change misfiring.
 *
 * The row counts reproduce DGG_getData's row block (:826-931):
 *   - `equalityRows`: `|rowUpper - rowLower| <= DGG_BOUND_THRESH` (1e-6). These are
 *     skipped outright by the slack loop at :1076 when mode == 0, i.e. for every
 *     tableau row.
 *   - `basicRows`: artificial basic. `activeSlackRows` = nonbasic and not equality,
 *     which is the set of rows that can contribute a slack coefficient to a
 *     tableau row -- an upper bound on the slack tail's length, not an exact count,
 *     because the BTRAN decides which of them are actually hit.
 *   - `negativeSlackRows`: `data->x[ncol+i] < -DGG_NULL_SLACK` (1e-5), the
 *     condition the generator itself flags as "row has negative slack" at :867.
 *     A correctly captured optimal basis should give zero; a nonzero count in a
 *     fixture means the loader reconstructed a different LP, so the bench treats
 *     it as a loud warning rather than a datum.
 *   - `integerSlackRows`: the `DGG_setIsInteger` test at :898-923 -- integral row
 *     bound(s) *and* every entry an integral coefficient on an integer column,
 *     all with `frac_part(...) > DGG_INTEGRALITY_THRESH` (1e-10) as the
 *     disqualifier. Note the asymmetry documented on cbcTwomirFixtureFracPart.
 *
 * `formBaseRows` / `formIntRows` size the formulation stage. It walks
 * `min(nrow, info.formulation_rows)` rows (:1617), builds each row plus its slack
 * into a base (:1117-1171), and `DGG_generateFormulationCutsFromBase` gives up at
 * :1676 unless at least one participant is integer (`tot_int == 0 goto CLEANUP`).
 * `formIntRows` counts the rows that pass that test, and it is what makes
 * `twomirCandidates == 0` an invalid skip on its own.
 */
static void cbcTwomirFixtureCount(const OsiSolverInterface *si,
  const CoinWarmStartBasis *ws, double threshold, int formulationRows,
  CbcTwomirFixtureCounts &c)
{
  memset(&c, 0, sizeof(c));

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

  // A basis written for a different model would index out of range; the caller
  // checks this too and skips, but refuse rather than count nonsense.
  if (ws->getNumStructural() != numberColumns || ws->getNumArtificial() != numberRows)
    return;

  const CoinPackedMatrix *byCol = si->getMatrixByCol();
  const int *columnLength = byCol ? byCol->getVectorLengths() : NULL;

  // The generator's own integer test is `si->getColType()[j] != 0`, which is
  // exactly `!isContinuous(j)`; isContinuous is used here for the same reason the
  // .ctype sidecar uses it -- it does not depend on the current bounds.
  std::vector< char > isInt(numberColumns, 0);

  for (int j = 0; j < numberColumns; ++j) {
    // CglTwomir.cpp:110, verbatim thresholds.
    if (colLower[j] < -1.0e20 && colUpper[j] > 1.0e20)
      ++c.freeColumns;
    if (colUpper[j] == colLower[j])
      ++c.fixedColumns;
    const bool integer = !si->isContinuous(j);
    isInt[j] = integer ? 1 : 0;
    if (integer)
      ++c.integerColumns;

    if (ws->getStructStatus(j) == CoinWarmStartBasis::basic) {
      if (integer) {
        ++c.basicIntegers;
        const double frac = cbcTwomirFixtureFracPart(colsol[j]);
        if (frac >= threshold && frac <= 1.0 - threshold)
          ++c.twomirCandidates;
      }
    } else {
      // No fixedness filter: CglTwomir.cpp:1057 has none. See the comment above.
      ++c.nonbasicStructural;
      if (columnLength)
        c.nonbasicStructuralNz += (double)columnLength[j];
    }
  }

  c.tabRowWork = (double)c.twomirCandidates * c.nonbasicStructuralNz;
  c.outerLoopWork = (double)c.twomirCandidates * (double)numberColumns;

  const CoinPackedMatrix *byRow = si->getMatrixByRow();
  if (!byRow)
    return;
  const int *column = byRow->getIndices();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *rowLength = byRow->getVectorLengths();
  const double *rowElements = byRow->getElements();

  const double infinity = si->getInfinity();
  c.formBaseRows = numberRows < formulationRows ? numberRows : formulationRows;
  if (c.formBaseRows < 0)
    c.formBaseRows = 0;

  for (int i = 0; i < numberRows; ++i) {
    const bool rowBasic = ws->getArtifStatus(i) == CoinWarmStartBasis::basic;
    // DGG_BOUND_THRESH, CglTwomir.cpp:834.
    const bool equality = fabs(rowUpper[i] - rowLower[i]) <= 1.0e-6;
    const bool above = rowUpper[i] < infinity;
    const bool below = rowLower[i] > -infinity;
    if (rowBasic)
      ++c.basicRows;
    if (equality)
      ++c.equalityRows;
    if (!rowBasic && !equality)
      ++c.activeSlackRows;

    // Row activity at colsol, exactly as CglTwomir.cpp:853-859 builds it, then the
    // slack value at :862-865.
    double activity = 0.0;
    for (CoinBigIndex k = rowStart[i]; k < rowStart[i] + rowLength[i]; ++k)
      activity += rowElements[k] * colsol[column[k]];
    const double slack = above ? rowUpper[i] - activity : activity - rowLower[i];
    if (slack < -1.0e-5) // DGG_NULL_SLACK
      ++c.negativeSlackRows;

    if (rowLength[i] > c.rowNzMax)
      c.rowNzMax = rowLength[i];

    // The integer-slack test, CglTwomir.cpp:898-923. DGG_INTEGRALITY_THRESH 1e-10.
    bool slackIsInteger = true;
    if (above) {
      if (cbcTwomirFixtureFracPart(rowUpper[i]) > 1.0e-10)
        slackIsInteger = false;
      else if (below && cbcTwomirFixtureFracPart(rowLower[i]) > 1.0e-10)
        slackIsInteger = false;
    } else if (cbcTwomirFixtureFracPart(rowLower[i]) > 1.0e-10) {
      slackIsInteger = false;
    }
    bool anyIntegerParticipant = false;
    for (CoinBigIndex k = rowStart[i]; k < rowStart[i] + rowLength[i]; ++k) {
      const int j = column[k];
      if (isInt[j])
        anyIntegerParticipant = true;
      if (slackIsInteger
        && (cbcTwomirFixtureFracPart(rowElements[k]) > 1.0e-10 || !isInt[j]))
        slackIsInteger = false;
    }
    if (slackIsInteger)
      ++c.integerSlackRows;

    // The formulation stage's own gate: DGG_getFormulaConstraint appends the slack
    // for a non-equality row (:1163-1170), and
    // DGG_generateFormulationCutsFromBase bails unless some participant is integer
    // (:1676).
    if (i < c.formBaseRows
      && (anyIntegerParticipant || (!equality && slackIsInteger)))
      ++c.formIntRows;
  }
}

/**
 * Write the provenance file.
 *
 * Three blocks: the generic ten (same names, order and formats as the Gomory dump
 * so a scan script can share a parser), the Twomir sizing counts, and the call
 * shape plus CBC's knob values. See the header comment for the derivation of every
 * `cbc*` and `eff*` key, and for why a bench that guesses them measures a call CBC
 * never makes.
 */
static bool cbcTwomirFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  int paddedColumns, int integerColumns, const CbcTwomirFixtureCounts &c, int pass,
  int options, int inTree, int level, int formulationRows, double awayAtRoot,
  double away, double cutoff, double maxSeconds, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;
  fprintf(fp, "tag %s\n", tag);
  fprintf(fp, "rows %d\n", si->getNumRows());
  fprintf(fp, "cols %d\n", si->getNumCols());
  fprintf(fp, "elements %d\n", si->getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si->isProvenOptimal());
  // The dumping run's budget, NOT a limit the Twomir call honoured: CbcCutGenerator
  // does call setMaxSeconds on this generator, but CglTwomir never reads it (no
  // maxSeconds / CoinCpuTime / CoinWallclockTime / getTimeOfDay / clock reference
  // anywhere in CglTwomir.cpp or .hpp), so the call has no internal deadline and a
  // replay is reproducible.
  fprintf(fp, "maxSeconds %.6f\n", maxSeconds);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index. `rows` above is
  // the captured count, so a consumer loading rows+1 deletes the last row.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());

  // What CglTwomir would look at: see cbcTwomirFixtureCount.
  // Nonzero freeColumns means generateCuts returns at CglTwomir.cpp:118 with no
  // cuts at all; the dump refuses to write such a fixture, so this is always 0
  // here and is recorded so the loader can assert it.
  fprintf(fp, "freeColumns %d\n", c.freeColumns);
  fprintf(fp, "fixedColumns %d\n", c.fixedColumns);
  fprintf(fp, "basicIntegers %d\n", c.basicIntegers);
  fprintf(fp, "twomirCandidates %d\n", c.twomirCandidates);
  fprintf(fp, "nonbasicStructural %d\n", c.nonbasicStructural);
  fprintf(fp, "nonbasicStructuralNz %.0f\n", c.nonbasicStructuralNz);
  // The leading term in tableau-stage cost: the dense loop at
  // CglTwomir.cpp:1056-1067 runs once per candidate over every element of every
  // nonbasic structural column. Precomputed so a ranking script need not multiply,
  // and as a double because the product overflows int on the larger instances.
  fprintf(fp, "tabRowWork %.0f\n", c.tabRowWork);
  // The floor the same loop cannot go below: the outer `for j < ncol` sweep must
  // stay dense and in increasing j to preserve emission order (CoinSort_2 at :1099
  // re-sorts only the slack tail). A fixture where this dominates tabRowWork
  // cannot be improved by sparsifying the inner loop.
  fprintf(fp, "outerLoopWork %.0f\n", c.outerLoopWork);
  fprintf(fp, "basicRows %d\n", c.basicRows);
  fprintf(fp, "equalityRows %d\n", c.equalityRows);
  fprintf(fp, "activeSlackRows %d\n", c.activeSlackRows);
  fprintf(fp, "negativeSlackRows %d\n", c.negativeSlackRows);
  fprintf(fp, "integerSlackRows %d\n", c.integerSlackRows);
  fprintf(fp, "formBaseRows %d\n", c.formBaseRows);
  fprintf(fp, "formIntRows %d\n", c.formIntRows);
  fprintf(fp, "rowNzMax %d\n", c.rowNzMax);

  // The call shape CBC used, so a bench can reproduce it rather than guess.
  fprintf(fp, "infoPass %d\n", pass);
  fprintf(fp, "infoOptions %d\n", options);
  fprintf(fp, "infoInTree %d\n", inTree);
  fprintf(fp, "infoLevel %d\n", level);
  fprintf(fp, "infoFormulationRows %d\n", formulationRows);
  fprintf(fp, "cutoff %.15g\n", cutoff);

  // CBC's Twomir configuration (CbcSolverCutSetup.cpp:359-418). It differs from
  // CglTwomir's constructor defaults on `away`/`awayAtRoot` in the *opposite*
  // direction from Gomory: the library default is 0.0005 (CglTwomir.cpp:570) and
  // CBC widens it to 0.005, a tenfold narrower candidate window. The values below
  // come from the `#ifdef MORE_CUTS` branch at :365-367, which is live only
  // because `#define MORE_CUTS` sits at :139 inside the *Gomory* block and is
  // never undefined.
  fprintf(fp, "cbcAwayAtRoot %.15g\n", awayAtRoot);
  fprintf(fp, "cbcAway %.15g\n", away);
  fprintf(fp, "cbcMaxElements 250\n");
  // The constructor's value, because setMaxElementsRoot is never called anywhere
  // in Cbc -- and it is max_elements_root_, not max_elements_, that :296 selects
  // at the root.
  fprintf(fp, "cbcMaxElementsRoot 50000\n");
  fprintf(fp, "cbcTwomirType 0\n");
  fprintf(fp, "cbcTMin 1\n");
  fprintf(fp, "cbcTMax 1\n");
  fprintf(fp, "cbcQMin 1\n");
  fprintf(fp, "cbcQMax 1\n");
  fprintf(fp, "cbcAMax 2\n");
  fprintf(fp, "cbcDoMir 1\n");
  fprintf(fp, "cbcDo2mir 1\n");
  fprintf(fp, "cbcDoTab 1\n");
  fprintf(fp, "cbcDoForm 1\n");
  fprintf(fp, "cbcHowOften -98\n");
  fprintf(fp, "cbcAccuracyFlag 4\n");
  fprintf(fp, "cbcSwitches 1\n");
  // switches[...] = 1 becomes setSwitchOffIfLessThan(1), and
  // CbcCutGenerator.cpp:1584-1592 then disables TwoMirCuts for the whole solve if
  // this very call returns fewer than 1 row cut. So the population this dump
  // defines is exactly the population the kill switch judges.
  fprintf(fp, "cbcSwitchOffIfLessThan 1\n");
  // CglTwomir's own generator, seeded at CglTwomir.cpp:569. info.randomNumberGenerator
  // is filled in by CbcCutGenerator but CglTwomir ignores it.
  fprintf(fp, "cbcSeed 987654321\n");

  // What the knobs ACTUALLY become inside the call, so the bench can assert it
  // reproduced the configuration instead of re-deriving it.
  // max_elements is overwritten with getNumCols() at :301-302 when
  // (!info.pass || (info.options & 32)) -- so at pass 0 it is neither 250 nor 50000.
  fprintf(fp, "effMaxElements %d\n",
    (!pass || (options & 32)) && !inTree ? si->getNumCols()
                                         : (inTree ? 250 : 50000));
  // data->gomory_threshold, :297. Note it is awayAtRoot_ outright at the root, not
  // min(away_, awayAtRoot_) as CglGomory does.
  fprintf(fp, "effGomoryThreshold %.15g\n", inTree ? away : awayAtRoot);
  // The tableau stage runs only when do_tab_ && info.level < 1 && info.pass < 6
  // (:314). Pass 6 onward is the formulation stage alone -- a structurally
  // different call, not a cheaper one.
  fprintf(fp, "effDoTabActive %d\n", (level < 1 && pass < 6) ? 1 : 0);
  fclose(fp);
  return true;
}

/**
 * Dump one fixture. Returns true if files were written.
 *
 * Six skip criteria. Four are shared with the Gomory dump; the last two are
 * Twomir's own and are the reason this population is not that one:
 *
 *   - `freeColumns > 0`: CglTwomir.cpp:109-118 counts columns free in both
 *     directions and returns outright if there are any, producing nothing. A
 *     fixture like that is not a slow Twomir case, it is a no-op.
 *   - `twomirCandidates == 0 && formIntRows == 0`, and only then. Unlike Gomory,
 *     no-candidates alone is NOT a valid skip: the formulation stage (:317-320)
 *     works from formulation rows and needs no fractional basic integer, so a
 *     Gomory-shaped gate would silently drop that whole population.
 */
static bool cbcDumpTwomirFixture(OsiSolverInterface *si, const char *tag = "twomir",
  int pass = 0, int options = 0, int inTree = 0, int formulationRows = -1,
  double maxSeconds = 0.0, int level = 0)
{
  if (!si)
    return false;

  const std::string name = cbcTwomirFixtureName(si);

  if (si->getNumCols() == 0 || si->getNumRows() == 0) {
    printf("[twomirfixture] %s.%s: SKIP (empty model: %d rows, %d cols)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols());
    fflush(stdout);
    return false;
  }

  if (!si->isProvenOptimal()) {
    printf("[twomirfixture] %s.%s: SKIP (LP not proven optimal; CglTwomir's "
           "needsOptimalBasis() is true and DGG_generateTabRowCuts factorizes "
           "this very basis, so a fixture without one is not replayable)\n",
      name.c_str(), tag);
    fflush(stdout);
    return false;
  }

  CoinWarmStart *rawWs = si->getWarmStart();
  const CoinWarmStartBasis *ws
    = dynamic_cast< const CoinWarmStartBasis * >(rawWs);
  if (!ws) {
    printf("[twomirfixture] %s.%s: SKIP (no CoinWarmStartBasis from the solver; "
           "DGG_getData does the same dynamic_cast at CglTwomir.cpp:711 and would "
           "dereference null)\n",
      name.c_str(), tag);
    fflush(stdout);
    delete rawWs;
    return false;
  }

  if (ws->getNumStructural() != si->getNumCols()
    || ws->getNumArtificial() != si->getNumRows()) {
    printf("[twomirfixture] %s.%s: SKIP (basis dimensions %d/%d do not match the "
           "model's %d cols / %d rows)\n",
      name.c_str(), tag, ws->getNumStructural(), ws->getNumArtificial(),
      si->getNumCols(), si->getNumRows());
    fflush(stdout);
    delete rawWs;
    return false;
  }

  // CBC's effective root threshold. CglTwomir.cpp:297 uses awayAtRoot_ outright
  // when !info.inTree -- no min() with away_, unlike CglGomory -- and CBC sets
  // awayAtRoot_ = 0.005, away_ = 0.01 (CbcSolverCutSetup.cpp:365-367).
  const double awayAtRoot = 0.005;
  const double away = 0.01;
  const double threshold = inTree ? away : awayAtRoot;

  const int effFormulationRows
    = formulationRows >= 0 ? formulationRows : si->getNumRows();

  CbcTwomirFixtureCounts c;
  cbcTwomirFixtureCount(si, ws, threshold, effFormulationRows, c);
  delete rawWs;

  if (c.freeColumns) {
    printf("[twomirfixture] %s.%s: SKIP (%d free column(s): CglTwomir.cpp:109-118 "
           "returns without generating anything when any column has "
           "lb < -1e20 && ub > 1e20, so this instance has no Twomir call to "
           "measure)\n",
      name.c_str(), tag, c.freeColumns);
    fflush(stdout);
    return false;
  }

  if (c.twomirCandidates == 0 && c.formIntRows == 0) {
    printf("[twomirfixture] %s.%s: SKIP (nothing for either stage; %d rows, "
           "%d cols, %d basic integer column(s) -- none has frac_part(colsol) "
           "inside (%g, %g) so the tableau stage is empty, and none of the %d "
           "formulation base row(s) has an integer participant so the formulation "
           "stage bails at CglTwomir.cpp:1676)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols(), c.basicIntegers,
      threshold, 1.0 - threshold, c.formBaseRows);
    fflush(stdout);
    return false;
  }

  double cutoff = COIN_DBL_MAX;
  si->getDblParam(OsiDualObjectiveLimit, cutoff);

  const std::string dir = cbcTwomirFixtureDir();
  cbcTwomirFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcTwomirFixtureWriteMps(si, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk
    = cbcTwomirFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const bool metaOk = cbcTwomirFixtureWriteMeta(si, tag, paddedColumns,
    integerColumns, c, pass, options, inTree, level, effFormulationRows,
    awayAtRoot, away, cutoff, maxSeconds, stem + ".meta");
  const bool basOk = cbcTwomirFixtureWriteBasis(si, stem + ".bas");
  const bool solOk = cbcTwomirFixtureWriteSol(si, stem + ".sol");

  printf("[twomirfixture] %s.%s: %s rows=%d cols=%d int=%d fixed=%d padded=%d "
         "cand=%d basicInt=%d nbStruct=%d nbStructNz=%.0f tabWork=%.0f "
         "outerWork=%.0f basicRows=%d eqRows=%d activeSlack=%d negSlack=%d "
         "intSlack=%d formBase=%d formInt=%d rowNzMax=%d obj=%.15g pass=%d "
         "options=%d (mps=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, c.fixedColumns,
    paddedColumns, c.twomirCandidates, c.basicIntegers, c.nonbasicStructural,
    c.nonbasicStructuralNz, c.tabRowWork, c.outerLoopWork, c.basicRows,
    c.equalityRows, c.activeSlackRows, c.negativeSlackRows, c.integerSlackRows,
    c.formBaseRows, c.formIntRows, c.rowNzMax, si->getObjValue(), pass, options,
    mpsRc, (int)ctOk, (int)metaOk, (int)basOk, (int)solOk);
  fflush(stdout);

  // The basis is required, not optional, so unlike the Probing dump it counts
  // toward success.
  return mpsRc == 0 && ctOk && metaOk && basOk;
}

#endif // CBC_DUMP_TWOMIR_FIXTURE
#endif // CbcTwomirFixtureDump_H
