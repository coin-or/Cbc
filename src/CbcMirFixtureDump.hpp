/**
 * Fixture dumping for CglMixedIntegerRounding2 experiments.
 *
 * @file CbcMirFixtureDump.hpp
 * @brief dump the exact inputs of a CglMixedIntegerRounding2 call, for offline replay
 *
 * Same motivation as CbcClqFixtureDump.hpp, CbcZeroHalfFixtureDump.hpp,
 * CbcProbingFixtureDump.hpp, CbcGomoryFixtureDump.hpp and
 * CbcTwomirFixtureDump.hpp: the call worth optimizing only happens deep inside a
 * solve, after preprocessing and the root LP, so experimenting on it costs a full
 * CBC run per measurement. This captures everything the call depends on:
 *
 *   <name>.mir.mps.gz   the problem as it stands at the call
 *   <name>.mir.bas      the LP basis -- required, though for a reason unlike Gomory's
 *   <name>.mir.sol      the LP solution, for cross-checking the warm start
 *   <name>.mir.ctype    which columns are integer (MPS cannot say: see below)
 *   <name>.mir.meta     key/value provenance, read back by the benchmark
 *
 * THE BASIS IS REQUIRED, BUT NOT BECAUSE ANYTHING FACTORIZES IT.
 * `CglMixedIntegerRounding2` does **not** override `needsOptimalBasis()` -- the
 * .hpp declares only `generateCuts`, `clone`, `refreshSolver` and `generateCpp` as
 * virtual overrides -- and it never touches `getBInvARow`, `getBInvRow` or a
 * `CoinFactorization`. It reads `getColSolution()`, `getRowActivity()`,
 * `getColUpper/Lower()`, `getRowSense()`, `getRightHandSide()` and
 * `getColType()`, and nothing else. So the Twomir argument ("the tableau row comes
 * out of a factorization of this basis") does not apply here at all.
 *
 * The basis is nonetheless the input, by a second and much less obvious route:
 * `mixIntRoundPreprocess` resolves every range row into a one-sided row **using
 * the current row activity**. CglMixedIntegerRounding2.cpp:502-518:
 *
 *     const double* rowActivity = si.getRowActivity();
 *     ...
 *     if (sense_[iRow]=='R') {
 *       if (rowActivity[iRow]-rowLower[iRow] < rowUpper[iRow]-rowActivity[iRow]) {
 *         RHS_[iRow]=rowLower[iRow];  sense_[iRow]='G';
 *       } else {
 *         RHS_[iRow]=rowUpper[iRow];  sense_[iRow]='L';   // ties go to 'L'
 *       }
 *     }
 *
 * `sense_` and `RHS_` then feed `determineRowType` (:520-524), so a different LP
 * point gives a different *row classification* -- different ROW_MIX / ROW_CONT /
 * ROW_INT sets, hence a different starting-row list, hence a different cut set.
 * A replay that starts from a cold LP is therefore not a weaker MIR2 fixture but a
 * different one, and this dump refuses to write a fixture whose LP is not proven
 * optimal. (An instance with no range rows would in principle replay from any
 * feasible point; `rangeRows` is recorded so that can be told apart, but the .bas
 * is written unconditionally because `getRowActivity()` is read unconditionally.)
 *
 * WHY THE PAYLOAD IS AGAIN JUST THE GENERIC FIVE. `refreshSolver()` is a no-op
 * body for this generator, and every array the workhorse touches is re-read from
 * the solver on each call. The one piece of state that *does* persist across calls
 * is the preprocessing result (`sense_`, `RHS_`, `rowTypes_`, `vubs_`, `vlbs_`,
 * `indRow*_`, `integerType_`, `numRows_`, `numCols_`, `doneInitPre_`) -- but CBC
 * calls `setDoPreproc(1)` (CbcSolverCutSetup.cpp:342), which makes
 * `generateCuts:64-67` recompute all of it on **every** call. So there is no
 * cross-call staleness to capture, and no analogue of the clique dump's
 * conflict-graph trap. MIR2 also draws no random numbers and reads no clock: there
 * is no `CoinCpuTime` / `CoinWallclockTime` / `maxSeconds` reference anywhere in
 * CglMixedIntegerRounding2.{cpp,hpp}, so a replay is reproducible.
 *
 * WHY THE GOMORY OR TWOMIR POPULATION CANNOT BE REUSED. The five *files* are
 * generator-agnostic and interchangeable; the *set of instances* is not.
 *
 *   - MIR2 HAS NO INSTANCE-WIDE DISQUALIFIER. Twomir returns outright when any
 *     column is free in both directions (CglTwomir.cpp:109-118), which is the
 *     single biggest reason its population differs from Gomory's. MIR2 has nothing
 *     like it, so this population should be *larger* than either.
 *   - ITS ONLY STRUCTURAL PRECONDITION IS `startingRows > 0`. The outer loop at
 *     CglMixedIntegerRounding2.cpp:866 runs `numRowMix_ + numRowContVB_ +
 *     numRowInt_` times, so when that sum is zero the generator does a full O(nnz)
 *     of preprocessing and matrix copying and then produces nothing. That is the
 *     one skip arm below that has no counterpart in the other dumps.
 *   - `rowsCont` ALONE IS NOT A CONTRIBUTOR. Only continuous rows that contain at
 *     least one column carrying a variable bound qualify (:700-715), so an
 *     instance can have hundreds of ROW_CONT rows and contribute none of them.
 *     Both counts are recorded; `rowsContVB` is the one in `startingRows`.
 *
 * The `.ctype` sidecar exists because MPS cannot express "fixed and integer": a
 * column with lb == ub takes writeMps's " FX " branch, which has no integer form,
 * and reads back continuous. For MIR2 a lost marker is worse than a weakening,
 * because integrality is what *classifies the rows*:
 *
 *   - `integerType_ = CoinCopyOfArray(si.getColType(), numCols_)` (:475), and
 *     `determineRowType` counts every entry as integer or continuous from it
 *     (:743-754). One lost marker can turn a ROW_INT into a ROW_MIX, a ROW_MIX
 *     into a ROW_CONT, or a two-variable ROW_VARUB into something unrecognised --
 *     changing the starting-row list, not merely one cut.
 *   - the variable-bound extraction at :654-664 picks the integer entry as `y` and
 *     the continuous entry as `x`; a lost marker swaps the roles and writes a
 *     nonsense `vubs_[x].setVal(-yCoef/xCoef)`.
 *   - `boundSubstitution` sends integer columns to the mixed knapsack and
 *     continuous ones through bound substitution (:1307-1311), so a lost marker
 *     moves a column from one side of the c-MIR to the other.
 *
 * WHAT THE .meta RECORDS ABOUT THE CALL, AND WHY EACH FIELD IS LOAD-BEARING. A
 * bench that calls generateCuts with a default-constructed
 * CglMixedIntegerRounding2 and a default CglTreeInfo measures a configuration CBC
 * never runs. For MIR2 the single most important fact is the one that is easiest to
 * get wrong:
 *
 *   - `cbcMaxAggr 1`. CBC constructs `CglMixedIntegerRounding2 mixedGen(1, true, 1)`
 *     (CbcSolverCutSetup.cpp:341) and only overrides MAXAGGR_ when
 *     `mixedRoundStrategy != 1` (:343-344) -- and `mixedRoundStrategy_` **defaults
 *     to 1** (Cbc/src/CbcSolver.hpp:206, re-asserted at CbcSolver.cpp:10926). So
 *     CBC runs MIR2 with `MAXAGGR_ == 1`, and the row aggregation the generator is
 *     named for **never happens**: the `iAggregate` loop at :879 runs exactly once
 *     with `iAggregate == 0` and takes the `copyRowSelected` arm at :880, so
 *     `selectRowToAggregate` (:1124) and `aggregateRow` (:1209) are never called.
 *     Any experiment aimed at the aggregation heuristic under CBC's defaults is an
 *     experiment on dead code.
 *   - `cbcMultiply 1` -> `effUpperLimit 2` (:839-842). Each starting row is
 *     separated twice, once as-is and once negated (:939-948). The two directions
 *     do not necessarily both reach `cMirSeparation`: the continuous-variables-in-s
 *     test inside `boundSubstitution` is sign-dependent (:1348, :1363), so `sStar`
 *     can be empty in one direction and not the other. `cmirCalls` below counts the
 *     directions that actually get through.
 *   - `cbcCriterion 1` -> `isLowerSubst` (:1238-1241) is the sign-independent
 *     `xlp - LB < UB - xlp`. Criteria 2 and 3 are implemented and never exercised
 *     by CBC; the bench can select them, which is why the sizing counts below take
 *     the criterion as a parameter rather than hardcoding 1.
 *   - `cbcDoPreproc 1` (:342). Preprocessing runs on every call, so the O(nnz)
 *     classification pass is part of the per-call cost, not a one-off.
 *   - `infoOptions 0`, re-derived for MIR2 rather than inherited from the other
 *     dumps. Every bit CbcCutGenerator.cpp:309-353 can set is clear at the first
 *     root pass: MIR2's howOften is `translate[CGIfMove] == -98`, which is neither
 *     `< -900` nor `< -1900`, so bits 8 and 16 stay clear; `strongCuts` is false on
 *     pass 1 because the cut_obj history is still -COIN_DBL_MAX (bit 32); fullScan
 *     is 2 at the root (bit 128); nothing in Cbc/src sets
 *     `moreSpecialOptions() & 16384` (bit 256); no parentModel (bit 512); not
 *     must-call-again (bit 1024). Bits 4 and 8 matter here for one specific reason:
 *     `generateCuts:215-219` marks the produced cuts globally valid exactly when
 *     `!info.inTree && ((options&4) || ((options&8) && !info.pass))`, and with
 *     options 0 that test is false -- so the root pass-0 cuts are *not* marked
 *     globally valid, and a bench that passed a nonzero options would produce cuts
 *     that differ in a field the comparison checks.
 *   - `cbcSwitches 0` (:347), because `lagrangeanFlag` is 0 by default. That
 *     becomes `setSwitchOffIfLessThan(0)`, a condition CbcCutGenerator's kill
 *     switch can never satisfy. **MIR2 has no first-pass kill switch**, unlike
 *     Twomir's `setSwitchOffIfLessThan(1)`. A scan that sees "MixedIntegerRounding2
 *     row present, Next run = disabled" is therefore looking at CbcCutGenerator's
 *     ordinary frequency tuning after the root, not at a first-call veto, and must
 *     record the string rather than attribute a mechanism to it.
 *   - `infoPass 0` is load-bearing in a way that has no Gomory or Twomir analogue.
 *     The wide-cut filter at :1001-1008 rejects any cut with
 *     `n > 0.8 * numCols_` -- but only `if (info_->pass || info_->inTree)`. At the
 *     root first pass it is **off**, so pass 0 accepts cuts that every later pass
 *     discards. This is MIR2's counterpart of CglGomory's `limitAtRoot_`, and it
 *     makes pass 0 both the widest and the most expensive call.
 *   - `maxSeconds` is the dumping run's budget, provenance only. CbcCutGenerator
 *     does call `setMaxSeconds(remaining)` on this generator, but the value lands on
 *     the base class and MIR2 never reads it -- no clock reference exists in the
 *     source -- so the call has no internal deadline.
 *
 * TWO OBSERVATIONS RECORDED HERE BECAUSE THIS IS WHERE THEY WERE FOUND. Neither is
 * fixed by a dump, and neither belongs in a performance change:
 *
 *   - `generateCuts:46-56` overwrites `MAXAGGR_` when it is negative and restores it
 *     from `saveMaxAggr` at :220 -- but the `MODIFY_LP==2` branch `return`s at :176,
 *     skipping the restore. So under `-mixedRoundStrategy -1` or `-2` the first root
 *     call leaves `MAXAGGR_` at 5 for the rest of the solve. Unreachable at CBC's
 *     default, since 1 is positive.
 *   - `boundSubstitution:1298-1305` handles `fabs(coefCol) < EPSILON_` by reading
 *     `colUpperBound[indCol]` / `colLowerBound[indCol]`, with no `indCol < numCols_`
 *     guard -- unlike the two branches around it. `copyRowSelected` gives every
 *     slack a coefficient of exactly +/-1, so the branch cannot fire for a slack at
 *     `MAXAGGR_ == 1`; at `MAXAGGR_ > 1` `aggregateRow` scales rows, and a scaled
 *     slack coefficient below 1e-6 would read one past the end of the bound arrays.
 *
 * Entirely behind CBC_DUMP_MIR_FIXTURE and off by default. Header-only static
 * functions, so no Makefile.am/Makefile.in changes are needed.
 *
 * Environment:
 *   CBC_MIR_FIXTURE_DIR    output directory
 *                          (default ~/instances/miplib/2017+spp/mirFixtures)
 *   CBC_MIR_FIXTURE_NAME   instance base name, when the MPS problem name is
 *                          unhelpful (drivers normally set this)
 **/

#ifndef CbcMirFixtureDump_H
#define CbcMirFixtureDump_H

#ifdef CBC_DUMP_MIR_FIXTURE

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

#include "CoinIndexedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

/**
 * Escape hatch for SKIP arm 5 (`startingRows == 0`). Set
 * `CBC_MIR_FIXTURE_ALLOW_EMPTY_START=1` to dump those instances instead of skipping
 * them; the resulting `.meta` carries `allowEmptyStart 1` so the sub-population is
 * self-identifying and never silently mixed into the main set.
 *
 * Why they are worth dumping rather than an oddity to discard. With zero starting
 * rows the generator still runs `mixIntRoundPreprocess` and still builds both O(nnz)
 * matrix copies, and only then executes the separation loop zero times -- so such a
 * fixture is a *pure fixed-cost probe*, in which 100% of the time is what OPT-A/B/C
 * remove and none of it is separation. Skipping them therefore biases the measured
 * speedup downward: in the first 38 results of the pass-1 population, 8 instances
 * skipped here and carried 13.54 s of root MIR2 time between them, `eilA101-2` alone
 * contributing 5.53 s across 50 passes for 0 cuts. They are kept separate because
 * they produce no cuts and so cannot participate in an exactness claim about cut
 * *content* -- only about cost.
 *
 * NOTE, and this cost time once: `startingRows == 0` is a pass-0 property and only a
 * pass-0 property. `numRows_ = si.getNumRows()` (CglMixedIntegerRounding2.cpp:459)
 * and `info.formulationRows` is never read, so from pass 1 on MIR2 classifies the cut
 * rows *other generators* added, and a single ROW_MIX cut row lifts `startingRows`
 * above zero. So a whole-solve cut count cannot falsify this arm: `nw04` skips here
 * and the generator table still credits it with a cut, yet at `-passC 1` it produces
 * 0 cuts in 0.037 s. To test a pass-0 property, re-run with a single cut pass.
 */
static bool cbcMirAllowEmptyStart()
{
  const char *e = getenv("CBC_MIR_FIXTURE_ALLOW_EMPTY_START");
  return e && *e && *e != '0';
}

/// mkdir -p, so a driver need not pre-create the tree.
static void cbcMirFixtureMkdirP(const std::string &path)
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
static std::string cbcMirFixtureDir()
{
  const char *env = getenv("CBC_MIR_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);

  const char *home = getenv("HOME");
  return std::string(home ? home : ".") + "/instances/miplib/2017+spp/mirFixtures";
}

/// Base name for this instance's files. The MPS problem name is the fallback, but
/// drivers should set CBC_MIR_FIXTURE_NAME: several instances share a problem
/// name, and some carry none at all.
static std::string cbcMirFixtureName(const OsiSolverInterface *si)
{
  const char *env = getenv("CBC_MIR_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);

  std::string probName;
  if (si->getStrParam(OsiProbName, probName) && !probName.empty())
    return probName;

  return std::string("instance");
}

/**
 * Write the MPS. formatType 2 is IEEE hex, which round-trips every double exactly;
 * plain MPS loses bits. For MIR2 that matters at three separate thresholds, all of
 * which are decisions rather than tolerances: `EPSILON_` (1e-6) decides in
 * `determineRowType` whether a coefficient is counted at all, and therefore what
 * *type* a row is; the same constant decides in `boundSubstitution` whether a
 * column joins the knapsack or is relaxed away (:1298); and the acceptance test at
 * :1027 rejects a cut outright when `largest > 1e8*smallest || largest > 1e7 ||
 * smallest < 1e-5`. A matrix perturbed in its last bits can move a row between
 * classes and so change the cut set structurally. numberAcross 1 for the same
 * reason.
 *
 * Empty columns are the one wrinkle: writeMps drops a column with no matrix entry,
 * which shifts every later column index and so invalidates the .ctype sidecar, the
 * .bas and the .sol. A single redundant final row covering those columns keeps them
 * at their original index; `paddedColumns` in the .meta says how many there were,
 * and the loader deletes the extra row.
 */
static int cbcMirFixtureWriteMps(const OsiSolverInterface *si,
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
    // would take the column indices with it. Note for MIR2 specifically that this
    // row is *ranged*, so `mixIntRoundPreprocess` will resolve it by activity like
    // any other range row; the loader deletes it before the generator ever sees it,
    // and `paddedColumns` in the .meta is what says the deletion is needed.
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
  // success*, leaving the replay to start cold. For MIR2 a cold start is not a
  // different vertex of the same problem, it is a different problem: the range-row
  // resolution at :506-517 reads getRowActivity(), so the row *types* change.
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
static bool cbcMirFixtureWriteBasis(OsiSolverInterface *si, const std::string &path)
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
static bool cbcMirFixtureWriteSol(const OsiSolverInterface *si, const std::string &path)
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
 * Binary/GeneralInteger from the current bounds. MIR2 itself uses
 * `si.getColType()` (:475) and only ever tests it for nonzero-ness, so the two
 * agree on the predicate that matters -- but see the header comment for the three
 * places a lost marker changes the row classification rather than merely weakening
 * a cut.
 */
static bool cbcMirFixtureWriteColTypes(const OsiSolverInterface *si,
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

/// The generator's RowType enum (CglMixedIntegerRounding2.hpp:96-119), reproduced
/// with the same ordinals so the .meta counts can be read against the source.
enum CbcMirFixtureRowType {
  CBC_MIR_ROW_UNDEFINED = 0,
  CBC_MIR_ROW_VARUB,
  CBC_MIR_ROW_VARLB,
  CBC_MIR_ROW_VAREQ,
  CBC_MIR_ROW_MIX,
  CBC_MIR_ROW_CONT,
  CBC_MIR_ROW_INT,
  CBC_MIR_ROW_OTHER
};

/**
 * `determineRowType` (CglMixedIntegerRounding2.cpp:722-810), reproduced rather than
 * approximated. Three details are easy to get wrong and each moves rows between
 * classes:
 *
 *   - The `+/- EPSILON_` band at :743-754 is a **dead zone, not a sign test**: a
 *     coefficient with `|a| <= EPSILON_` increments neither the integer nor the
 *     continuous count. So a row of nothing but tiny coefficients has
 *     `numInt == 0` and falls through to ROW_CONT at :796, not to ROW_OTHER.
 *   - The variable-bound test at :774 requires `|rhs| <= EPSILON_` as well as
 *     exactly one integer and one continuous entry, and an 'E' row that passes it
 *     becomes ROW_VAREQ while an 'E' row that fails it becomes ROW_MIX.
 *   - `numInt == 0` is tested (:796) **before** `numCon == 0` (:800), so the
 *     no-variables-at-all case is ROW_CONT.
 *
 * `isInt[]` is indexed by column and is `!isContinuous(j)`, which is the same
 * predicate as the generator's `integerType_[j] != 0`.
 */
static int cbcMirFixtureRowType(int rowLen, const int *ind, const double *coef,
  char sense, double rhs, const char *isInt, double epsilon)
{
  if (rowLen == 0 || fabs(rhs) > 1.0e20 || sense == 'N')
    return CBC_MIR_ROW_UNDEFINED;

  int numPosInt = 0, numNegInt = 0, numPosCon = 0, numNegCon = 0;
  for (int i = 0; i < rowLen; ++i) {
    if (coef[i] < -epsilon) {
      if (isInt[ind[i]])
        ++numNegInt;
      else
        ++numNegCon;
    } else if (coef[i] > epsilon) {
      if (isInt[ind[i]])
        ++numPosInt;
      else
        ++numPosCon;
    }
  }
  const int numInt = numNegInt + numPosInt;
  const int numCon = numNegCon + numPosCon;

  int rowType = CBC_MIR_ROW_UNDEFINED;
  if (numInt > 0 && numCon > 0) {
    if (numInt == 1 && numCon == 1 && fabs(rhs) <= epsilon) {
      switch (sense) {
      case 'L':
        rowType = numPosCon == 1 ? CBC_MIR_ROW_VARUB : CBC_MIR_ROW_VARLB;
        break;
      case 'G':
        rowType = numPosCon == 1 ? CBC_MIR_ROW_VARLB : CBC_MIR_ROW_VARUB;
        break;
      case 'E':
        rowType = CBC_MIR_ROW_VAREQ;
        break;
      default:
        break;
      }
    } else {
      rowType = CBC_MIR_ROW_MIX;
    }
  } else if (numInt == 0) {
    rowType = CBC_MIR_ROW_CONT;
  } else if (numCon == 0 && (sense == 'L' || sense == 'G')) {
    rowType = CBC_MIR_ROW_INT;
  } else {
    rowType = CBC_MIR_ROW_OTHER;
  }
  return rowType;
}

/// `isLowerSubst` (CglMixedIntegerRounding2.cpp:1231-1254), reproduced. Criterion 1
/// -- CBC's -- ignores the coefficient entirely, so the substitution branch is the
/// same for a row and its negation; criteria 2 and 3 test the sign, so under those
/// the two multiply directions can send *different* integer columns into the
/// knapsack. That is the reason this is a parameter here rather than a constant.
static bool cbcMirFixtureIsLowerSubst(int criterion, double infinity, double aj,
  double xlp, double LB, double UB)
{
  if (criterion == 1)
    return xlp - LB < UB - xlp;
  if (UB == infinity || xlp == LB)
    return true;
  if (LB == -infinity || xlp == UB)
    return false;
  if (criterion == 2)
    return aj < 0;
  return aj > 0;
}

/// Everything cbcMirFixtureCount produces. A struct rather than two dozen out
/// parameters, and the field names are exactly the .meta keys.
struct CbcMirFixtureCounts {
  int integerColumns;
  int fixedColumns;
  int rangeRows;
  int geRowsToFlip;
  int rowsUndefined;
  int rowsVarUB;
  int rowsVarLB;
  int rowsVarEQ;
  int rowsMix;
  int rowsCont;
  int rowsContVB;
  int rowsInt;
  int rowsOther;
  int startingRows;
  int vubCount;
  int vlbCount;
  int cmirCalls;
  int knapsackIntMax;
  double knapsackIntSum;
  double cmirWork;
  double nnzWork;
  int rowNzMax;
};

/**
 * Count what CglMixedIntegerRounding2 would actually work on, mirroring
 * `mixIntRoundPreprocess` (:453-717), `determineRowType` (:722-810),
 * `copyRowSelected` (:1085-1119) and `boundSubstitution` (:1261-1413) rather than
 * approximating them. Every count below is the generator's own notion of work.
 *
 * THE FIVE ROW-TYPE COUNTS ARE COMPUTED HERE, NOT MIRRORED FROM THE GENERATOR.
 * `mixIntRoundPreprocess` keeps only three counters unconditionally -- `numMIX`,
 * `numCONT`, `numINT` (:494-496). `numUNDEFINED`, `numVARUB`, `numVARLB`,
 * `numVAREQ` and `numOTHER` exist only under `#if CGL_DEBUG`, and CGL_DEBUG
 * defaults to 0 (CglMixedIntegerRounding2.hpp:26-28), so in a production build they
 * are not compiled at all. `rowsVarUB` / `rowsVarLB` / `rowsVarEQ` / `rowsOther` /
 * `rowsUndefined` below are therefore this dump's own tally of the same
 * classification, not a readback of a live counter.
 *
 * `rangeRows` and `geRowsToFlip` are the two solution- and shape-dependent inputs:
 *
 *   - `rangeRows` counts `getRowSense() == 'R'`. Each one is resolved to 'G' or 'L'
 *     by comparing the row activity to the two bounds (:508-517), **ties going to
 *     'L'**, and that resolved sense is what `determineRowType` classifies. This is
 *     the entire reason the fixture needs a .bas.
 *   - `geRowsToFlip` counts `rowUpper >= 1e50`, which is the `MODIFY_LP==2` swap
 *     test at :137-142. When it is zero, `nChanged == 0`, the deep copy is thrown
 *     away at :177-178 and the generator does **not** recurse -- so the call the
 *     bench reproduces at `--level=0` is structurally different from the one it
 *     reproduces on an instance with `>= 1` rows. A bench measuring the recursion
 *     overhead on a `geRowsToFlip == 0` fixture is measuring nothing.
 *
 * `startingRows = rowsMix + rowsContVB + rowsInt` is the outer loop bound at :854
 * and :866, and it is the generator's one structural precondition: at zero the loop
 * body never runs and MIR2 produces nothing after a full O(nnz) of preprocessing
 * and matrix copying. Note that `rowsCont` is **not** a term -- only the ROW_CONT
 * rows containing a column that carries a variable bound qualify (:700-715).
 *
 * `cmirCalls`, `knapsackIntSum`, `knapsackIntMax` and `cmirWork` come from
 * simulating each starting row through `copyRowSelected` and `boundSubstitution`,
 * for each of the `upperLimit` multiply directions:
 *
 *   - `copyRowSelected` copies the row's entries **unmodified** -- there is no sign
 *     normalisation -- and appends one slack entry at column index `numCols_` with
 *     coefficient +1 for an 'L' row and -1 for a 'G' row (:1110-1117). An 'E' row
 *     gets no slack. The `i == 1` direction then negates everything (:946).
 *   - `boundSubstitution` sends a column to the mixed knapsack when it is an
 *     unfixed integer model column (:1307-1311) or when it is a continuous column
 *     whose chosen bound is a *variable* bound, in which case the bound's integer
 *     variable enters instead (:1342, :1357). `numInt` is the resulting
 *     `CoinIndexedVector`'s element count, and an index enters that vector exactly
 *     once, the first time a value of magnitude at least COIN_INDEXED_TINY_ELEMENT
 *     lands on it (CoinIndexedVector::add) -- which is what the marker set below
 *     reproduces.
 *   - A direction reaches `cMirSeparation` only if `boundSubstitution` returns
 *     true, and it can fail four ways: a continuous column with both bounds
 *     infinite (:1328-1335), no continuous variable in s (:1388), an empty knapsack
 *     (:1398), or a knapsack integer column with a nonzero lower bound (:1407).
 *     The second of those is **sign-dependent** -- the lower-substitution branch
 *     admits a column to s only when `coefCol < -EPSILON_` and the upper branch
 *     only when `coefCol > EPSILON_` -- so one multiply direction can separate
 *     while the other does not. `cmirCalls` counts the directions that get through,
 *     which is the exact number of `cMirSeparation` entries per call.
 *
 * `cmirWork` is `sum over separating directions of numInt^2`, which is the leading
 * term in separation cost: the delta scan at :1487-1512 makes up to `numInt` calls
 * to `cMirInequality` and each is O(numInt), the delta/2,/4,/8 block at :1523-1541
 * adds three more, and the complemented-set loop at :1550-1578 adds up to
 * `complTSize <= numInt` more. Rank the slow tail by this.
 *
 * `nnzWork` is just the element count, recorded as the per-call O(nnz) **floor**
 * that the two matrix copies at :191-195 sit on regardless of how little
 * separation happens. It is the counterpart of Twomir's `outerLoopWork`: a fixture
 * where `nnzWork` dominates `cmirWork` cannot be improved by touching the
 * separation heuristic, and without publishing both such a fixture looks like a
 * change misfiring.
 *
 * `cmirWork` and `knapsackIntSum` are doubles because the products overflow int on
 * the larger instances.
 */
static void cbcMirFixtureCount(const OsiSolverInterface *si,
  const CoinWarmStartBasis *ws, int criterion, double epsilon, int undefined,
  int upperLimit, CbcMirFixtureCounts &c)
{
  memset(&c, 0, sizeof(c));

  const int numberColumns = si->getNumCols();
  const int numberRows = si->getNumRows();
  if (numberColumns == 0 || numberRows == 0)
    return;

  const double *colLower = si->getColLower();
  const double *colUpper = si->getColUpper();
  const double *rowLower = si->getRowLower();
  const double *rowUpper = si->getRowUpper();
  const double *colsol = si->getColSolution();
  const double *rowActivity = si->getRowActivity();
  const char *rowSense = si->getRowSense();
  const double *rhsIn = si->getRightHandSide();
  if (!colsol || !rowActivity || !rowSense || !rhsIn)
    return;

  // A basis written for a different model would make the .bas unusable on replay;
  // the caller checks this too and skips. MIR2 itself never reads the basis
  // statuses -- see the header comment on why the .bas is still the input.
  if (ws
    && (ws->getNumStructural() != numberColumns
      || ws->getNumArtificial() != numberRows))
    return;

  const CoinPackedMatrix *byRow = si->getMatrixByRow();
  if (!byRow)
    return;
  const double *rowElements = byRow->getElements();
  const int *column = byRow->getIndices();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *rowLength = byRow->getVectorLengths();

  c.nnzWork = (double)si->getNumElements();

  std::vector< char > isInt(numberColumns, 0);
  for (int j = 0; j < numberColumns; ++j) {
    if (!si->isContinuous(j)) {
      isInt[j] = 1;
      ++c.integerColumns;
    }
    if (colUpper[j] == colLower[j])
      ++c.fixedColumns;
  }

  // ---- mixIntRoundPreprocess, first pass: resolve ranges, classify rows.
  std::vector< char > sense(numberRows, 'N');
  std::vector< double > rhs(numberRows, 0.0);
  std::vector< char > rowType(numberRows, (char)CBC_MIR_ROW_UNDEFINED);

  for (int i = 0; i < numberRows; ++i) {
    sense[i] = rowSense[i];
    rhs[i] = rhsIn[i];
    if (rowLength[i] > c.rowNzMax)
      c.rowNzMax = rowLength[i];
    // The MODIFY_LP==2 swap test, :137-142, on the same bounds the copy would see.
    if (rowUpper[i] >= 1.0e50)
      ++c.geRowsToFlip;

    if (sense[i] == 'R') {
      ++c.rangeRows;
      // :508-517, ties to 'L'.
      if (rowActivity[i] - rowLower[i] < rowUpper[i] - rowActivity[i]) {
        rhs[i] = rowLower[i];
        sense[i] = 'G';
      } else {
        rhs[i] = rowUpper[i];
        sense[i] = 'L';
      }
    }

    const int t = cbcMirFixtureRowType(rowLength[i], column + rowStart[i],
      rowElements + rowStart[i], sense[i], rhs[i], isInt.empty() ? NULL : &isInt[0],
      epsilon);
    rowType[i] = (char)t;
    switch (t) {
    case CBC_MIR_ROW_UNDEFINED: ++c.rowsUndefined; break;
    case CBC_MIR_ROW_VARUB:     ++c.rowsVarUB;     break;
    case CBC_MIR_ROW_VARLB:     ++c.rowsVarLB;     break;
    case CBC_MIR_ROW_VAREQ:     ++c.rowsVarEQ;     break;
    case CBC_MIR_ROW_MIX:       ++c.rowsMix;       break;
    case CBC_MIR_ROW_CONT:      ++c.rowsCont;      break;
    case CBC_MIR_ROW_INT:       ++c.rowsInt;       break;
    default:                    ++c.rowsOther;     break;
    }
  }

  // ---- mixIntRoundPreprocess, second pass: extract the variable bounds, :622-690.
  // `xInd`/`yInd` start at 0 exactly as the generator's do (:651), so a row whose
  // every coefficient is within the epsilon band writes vubs_[0]/vlbs_[0] with a
  // 0/0 ratio. Reproduced rather than guarded, because it changes which columns
  // report a variable bound and therefore which ROW_CONT rows qualify below.
  std::vector< int > vubVar(numberColumns, undefined), vlbVar(numberColumns, undefined);
  std::vector< double > vubVal(numberColumns, 0.0), vlbVal(numberColumns, 0.0);

  for (int i = 0; i < numberRows; ++i) {
    const int t = rowType[i];
    if (t != CBC_MIR_ROW_VARUB && t != CBC_MIR_ROW_VARLB && t != CBC_MIR_ROW_VAREQ)
      continue;
    int xInd = 0, yInd = 0;
    double xCoef = 0.0, yCoef = 0.0;
    const CoinBigIndex startPos = rowStart[i];
    const CoinBigIndex stopPos = startPos + rowLength[i];
    for (CoinBigIndex k = startPos; k < stopPos; ++k) {
      if (fabs(rowElements[k]) > epsilon) {
        if (isInt[column[k]]) {
          yInd = column[k];
          yCoef = rowElements[k];
        } else {
          xInd = column[k];
          xCoef = rowElements[k];
        }
      }
    }
    const double val = -yCoef / xCoef;
    if (t == CBC_MIR_ROW_VARUB || t == CBC_MIR_ROW_VAREQ) {
      vubVar[xInd] = yInd;
      vubVal[xInd] = val;
    }
    if (t == CBC_MIR_ROW_VARLB || t == CBC_MIR_ROW_VAREQ) {
      vlbVar[xInd] = yInd;
      vlbVal[xInd] = val;
    }
  }

  for (int j = 0; j < numberColumns; ++j) {
    if (vubVar[j] != undefined)
      ++c.vubCount;
    if (vlbVar[j] != undefined)
      ++c.vlbCount;
  }

  // ---- The starting-row list, in the generator's own order (:883-891): the
  // ROW_MIX rows, then the ROW_CONT rows that carry a variable bound (:700-715),
  // then the ROW_INT rows.
  std::vector< int > starting;
  starting.reserve(c.rowsMix + c.rowsCont + c.rowsInt);
  for (int i = 0; i < numberRows; ++i) {
    if (rowType[i] == CBC_MIR_ROW_MIX)
      starting.push_back(i);
  }
  for (int i = 0; i < numberRows; ++i) {
    if (rowType[i] != CBC_MIR_ROW_CONT)
      continue;
    const CoinBigIndex jStart = rowStart[i];
    const CoinBigIndex jStop = jStart + rowLength[i];
    for (CoinBigIndex k = jStart; k < jStop; ++k) {
      const int indCol = column[k];
      if (vlbVar[indCol] != undefined || vubVar[indCol] != undefined) {
        starting.push_back(i);
        ++c.rowsContVB;
        break;
      }
    }
  }
  for (int i = 0; i < numberRows; ++i) {
    if (rowType[i] == CBC_MIR_ROW_INT)
      starting.push_back(i);
  }
  c.startingRows = c.rowsMix + c.rowsContVB + c.rowsInt;

  // ---- boundSubstitution, simulated per starting row per multiply direction.
  // The knapsack is a set of column indices; `mark` plus `touched` reproduces
  // CoinIndexedVector's "an index is appended the first time a value of magnitude
  // at least COIN_INDEXED_TINY_ELEMENT lands on it, and is never removed".
  const double infinity = si->getInfinity();
  const double tiny = COIN_INDEXED_TINY_ELEMENT;
  std::vector< char > mark(numberColumns, 0);
  std::vector< int > touched;
  touched.reserve(numberColumns);

  for (size_t s = 0; s < starting.size(); ++s) {
    const int iRow = starting[s];
    const CoinBigIndex jStart = rowStart[iRow];
    const CoinBigIndex jStop = jStart + rowLength[iRow];
    // The slack copyRowSelected would append: +1 for 'L', -1 for 'G', none for 'E'.
    double slackCoef = 0.0;
    if (sense[iRow] == 'L')
      slackCoef = 1.0;
    else if (sense[iRow] == 'G')
      slackCoef = -1.0;

    for (int dir = 0; dir < upperLimit; ++dir) {
      const double sign = dir == 0 ? 1.0 : -1.0;
      int numCont = 0;
      bool bailed = false;

      for (size_t t = 0; t < touched.size(); ++t)
        mark[touched[t]] = 0;
      touched.clear();

      for (CoinBigIndex k = jStart; k < jStop && !bailed; ++k) {
        const int indCol = column[k];
        const double coefCol = sign * rowElements[k];

        if (colLower[indCol] == colUpper[indCol])
          continue;                        // :1292-1296, fixed column removed
        if (fabs(coefCol) < epsilon)
          continue;                        // :1298-1305, relaxed away
        if (isInt[indCol]) {               // :1307-1311
          if (fabs(coefCol) >= tiny && !mark[indCol]) {
            mark[indCol] = 1;
            touched.push_back(indCol);
          }
          continue;
        }

        // Continuous model column: pick a bound, :1315-1368.
        const double LB = vlbVar[indCol] != undefined
          ? vlbVal[indCol] * colsol[vlbVar[indCol]]
          : colLower[indCol];
        const double UB = vubVar[indCol] != undefined
          ? vubVal[indCol] * colsol[vubVar[indCol]]
          : colUpper[indCol];
        if (LB == -1.0 * infinity && UB == infinity) {
          bailed = true;                   // :1328-1335, no mixed knapsack
          break;
        }
        if (cbcMirFixtureIsLowerSubst(criterion, infinity, coefCol, colsol[indCol],
              LB, UB)) {
          if (vlbVar[indCol] != undefined) {
            const int indVLB = vlbVar[indCol];
            const double v = coefCol * vlbVal[indCol];
            if (fabs(v) >= tiny && !mark[indVLB]) {
              mark[indVLB] = 1;
              touched.push_back(indVLB);
            }
          }
          if (coefCol < -epsilon)
            ++numCont;                     // :1348-1352
        } else {
          if (vubVar[indCol] != undefined) {
            const int indVUB = vubVar[indCol];
            const double v = coefCol * vubVal[indCol];
            if (fabs(v) >= tiny && !mark[indVUB]) {
              mark[indVUB] = 1;
              touched.push_back(indVUB);
            }
          }
          if (coefCol > epsilon)
            ++numCont;                     // :1363-1367
        }
      }

      // The slack, which copyRowSelected placed at index numCols_. It is never a
      // knapsack entry (the `indCol < numCols_` guards at :1292 and :1307 exclude
      // it); it can only join s, and only in one of the two directions, :1374-1378.
      if (!bailed && slackCoef != 0.0 && sign * slackCoef < -epsilon)
        ++numCont;

      if (bailed || numCont == 0)
        continue;                          // :1388
      const int numInt = (int)touched.size();
      if (numInt == 0)
        continue;                          // :1398
      bool zeroLower = true;
      for (int t = 0; t < numInt; ++t) {
        const int indCol = touched[t];
        if (fabs(colLower[indCol]) > epsilon) {
          zeroLower = false;               // :1402-1408
          break;
        }
      }
      if (!zeroLower)
        continue;

      ++c.cmirCalls;
      c.knapsackIntSum += (double)numInt;
      if (numInt > c.knapsackIntMax)
        c.knapsackIntMax = numInt;
      c.cmirWork += (double)numInt * (double)numInt;
    }
  }
}

/**
 * Write the provenance file.
 *
 * Three blocks: the generic ten (same names, order and formats as the Gomory and
 * Twomir dumps so a scan script can share a parser), the MIR2 sizing counts, and
 * the call shape plus CBC's knob values. See the header comment for the derivation
 * of every `cbc*` and `eff*` key, and for why a bench that guesses them measures a
 * call CBC never makes.
 */
static bool cbcMirFixtureWriteMeta(const OsiSolverInterface *si, const char *tag,
  int paddedColumns, int integerColumns, const CbcMirFixtureCounts &c, int pass,
  int options, int inTree, int level, int formulationRows, double cutoff,
  double maxSeconds, const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;
  fprintf(fp, "tag %s\n", tag);
  fprintf(fp, "rows %d\n", si->getNumRows());
  fprintf(fp, "cols %d\n", si->getNumCols());
  fprintf(fp, "elements %d\n", si->getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si->isProvenOptimal());
  // The dumping run's budget, NOT a limit the MIR2 call honoured: CbcCutGenerator
  // does call setMaxSeconds on this generator, but CglMixedIntegerRounding2 never
  // reads it -- there is no maxSeconds / CoinCpuTime / CoinWallclockTime /
  // getTimeOfDay / clock reference anywhere in its .cpp or .hpp -- so the call has
  // no internal deadline and a replay is reproducible.
  fprintf(fp, "maxSeconds %.6f\n", maxSeconds);
  // Nonzero means the .mps carries one extra, redundant final row so that this
  // many empty columns survive the write at their original index. `rows` above is
  // the captured count, so a consumer loading rows+1 deletes the last row.
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);
  fprintf(fp, "objSense %g\n", si->getObjSense());
  fprintf(fp, "objValue %.15g\n", si->getObjValue());

  // What CglMixedIntegerRounding2 would look at: see cbcMirFixtureCount.
  fprintf(fp, "fixedColumns %d\n", c.fixedColumns);
  // Range rows are the reason this fixture needs a .bas: each is resolved to 'G' or
  // 'L' from the row activity (:508-517, ties to 'L') before classification, so a
  // different LP point gives a different set of row types.
  fprintf(fp, "rangeRows %d\n", c.rangeRows);
  // The MODIFY_LP==2 swap count (:137-142). Zero means nChanged == 0, the deep copy
  // is discarded at :177-178 and the generator does NOT recurse -- a structurally
  // different call from the usual one, and one on which the recursion overhead
  // cannot be measured.
  fprintf(fp, "geRowsToFlip %d\n", c.geRowsToFlip);
  // The row classification. These five are computed by this dump, not read back
  // from the generator: only numMIX/numCONT/numINT are counted outside
  // `#if CGL_DEBUG`, and CGL_DEBUG is 0 by default.
  fprintf(fp, "rowsUndefined %d\n", c.rowsUndefined);
  fprintf(fp, "rowsVarUB %d\n", c.rowsVarUB);
  fprintf(fp, "rowsVarLB %d\n", c.rowsVarLB);
  fprintf(fp, "rowsVarEQ %d\n", c.rowsVarEQ);
  fprintf(fp, "rowsOther %d\n", c.rowsOther);
  fprintf(fp, "rowsMix %d\n", c.rowsMix);
  fprintf(fp, "rowsCont %d\n", c.rowsCont);
  // Only the ROW_CONT rows that contain a column carrying a variable bound
  // (:700-715). This, not rowsCont, is the term in startingRows.
  fprintf(fp, "rowsContVB %d\n", c.rowsContVB);
  fprintf(fp, "rowsInt %d\n", c.rowsInt);
  fprintf(fp, "vubCount %d\n", c.vubCount);
  fprintf(fp, "vlbCount %d\n", c.vlbCount);
  // The outer loop bound at :854/:866, and the generator's only structural
  // precondition: at zero the loop body never runs.
  fprintf(fp, "startingRows %d\n", c.startingRows);
  // Whether arm 5 was disabled for this dump. A fixture with startingRows == 0 is a
  // pure fixed-cost probe (preprocess + both O(nnz) copies, then a zero-iteration
  // separation loop) and belongs in its own labelled sub-population, not mixed in.
  fprintf(fp, "allowEmptyStart %d\n", cbcMirAllowEmptyStart() ? 1 : 0);
  // How many (starting row, multiply direction) pairs actually reach
  // cMirSeparation, i.e. how many times boundSubstitution returned true. The two
  // directions need not both get through: the continuous-variables-in-s test is
  // sign-dependent (:1348, :1363).
  fprintf(fp, "cmirCalls %d\n", c.cmirCalls);
  fprintf(fp, "knapsackIntMax %d\n", c.knapsackIntMax);
  fprintf(fp, "knapsackIntSum %.0f\n", c.knapsackIntSum);
  // The leading term in separation cost: sum of numInt^2 over separating
  // directions. cMirSeparation makes up to numInt calls to cMirInequality in the
  // delta scan (:1487-1512), three more at delta/2,/4,/8 (:1523-1541) and up to
  // complTSize <= numInt more in the complemented-set loop (:1550-1578), each O(numInt).
  // Rank the slow tail by this. A double because the product overflows int.
  fprintf(fp, "cmirWork %.0f\n", c.cmirWork);
  // The per-call O(nnz) floor the two matrix copies at :191-195 sit on regardless
  // of how little separation happens. A fixture where this dominates cmirWork
  // cannot be improved by touching the separation heuristic.
  fprintf(fp, "nnzWork %.0f\n", c.nnzWork);
  fprintf(fp, "rowNzMax %d\n", c.rowNzMax);

  // The call shape CBC used, so a bench can reproduce it rather than guess.
  fprintf(fp, "infoPass %d\n", pass);
  fprintf(fp, "infoOptions %d\n", options);
  fprintf(fp, "infoInTree %d\n", inTree);
  fprintf(fp, "infoLevel %d\n", level);
  // CbcCutGenerator fills this in from model_->numberRowsAtContinuous(); MIR2 never
  // reads info.formulation_rows. Recorded for parser compatibility with the Gomory
  // and Twomir .meta files, and as provenance.
  fprintf(fp, "infoFormulationRows %d\n", formulationRows);
  fprintf(fp, "cutoff %.15g\n", cutoff);

  // CBC's MIR2 configuration (CbcSolverCutSetup.cpp:338-348). For this generator
  // CBC's values happen to equal the alternate constructor's (1, true, 1) plus
  // setDoPreproc(1) -- but only because mixedRoundStrategy_ defaults to 1, which is
  // what makes the aggregation dead code. See the header comment.
  fprintf(fp, "cbcMaxAggr 1\n");
  fprintf(fp, "cbcMultiply 1\n");
  fprintf(fp, "cbcCriterion 1\n");
  fprintf(fp, "cbcDoPreproc 1\n");
  fprintf(fp, "cbcHowOften -98\n");
  fprintf(fp, "cbcAccuracyFlag 2\n");
  // switches[...] = 0 | (ALL_LAGRANGEAN * lagrangeanFlag), and lagrangeanFlag is 0
  // by default, so this becomes setSwitchOffIfLessThan(0) -- a condition
  // CbcCutGenerator's kill switch can never satisfy. MIR2 has NO first-pass kill
  // switch, unlike Twomir's 1. A "Next run = disabled" in the generator table is
  // ordinary frequency tuning after the root, not a first-call veto.
  fprintf(fp, "cbcSwitches 0\n");
  fprintf(fp, "cbcSwitchOffIfLessThan 0\n");
  // Cut acceptance: insertIfNotDuplicateAndClean(cut, 31, CoinAbsFltEq(1.0e-4)) at
  // :1044 and :865. Row cuts only -- there is no OsiColCut anywhere in the
  // generator -- and the near-duplicate test makes the *generation order* part of
  // the contract, since which of two similar cuts survives depends on which
  // arrived first.
  fprintf(fp, "cbcDuplicateTol 0.0001\n");

  // What the knobs ACTUALLY become inside the call, so the bench can assert it
  // reproduced the configuration instead of re-deriving it.
  // MULTIPLY_ true -> upperLimit 2 (:838-842): each starting row is separated twice,
  // once as-is and once negated.
  fprintf(fp, "effUpperLimit 2\n");
  // The wide-cut filter at :1001-1008 rejects cuts with n > 0.8*numCols_, but only
  // `if (info_->pass || info_->inTree)`. At the root first pass it is OFF, so this
  // call accepts cuts every later pass would discard.
  fprintf(fp, "effWideCutFilterOff %d\n", (!pass && !inTree) ? 1 : 0);
  // MAXAGGR_ == 1 leaves no reachable reader for coefByCol / rowInds / colStarts:
  // selectRowToAggregate (:904) and the SAFE_ROWS row-type clobber (:925-928) are
  // both in the `else` arm of the `iAggregate == 0` split at :880/:897, and the
  // SAFE_ROWS restore loop (:1053) starts at jAggregate = 1 and so runs zero times.
  // The matrixByCol transpose at :195 is therefore built and never read.
  fprintf(fp, "effMatrixByColDead 1\n");
  fclose(fp);
  return true;
}

/**
 * Dump one fixture. Returns true if files were written.
 *
 * Five skip criteria. The first four are shared with the Gomory and Twomir dumps;
 * the fifth is MIR2's own and is the only structural precondition it has:
 *
 *   - `startingRows == 0`: the outer loop at CglMixedIntegerRounding2.cpp:866 runs
 *     `numRowMix_ + numRowContVB_ + numRowInt_` times, so at zero the generator
 *     does a full O(nnz) of preprocessing and matrix copying and then produces
 *     nothing. There is no `freeColumns`-style instance-wide disqualifier here as
 *     there is in Twomir, and no `gomoryCandidates`-style candidate count either.
 *     Arm 5 alone has an escape hatch, `CBC_MIR_FIXTURE_ALLOW_EMPTY_START=1` -- see
 *     `cbcMirAllowEmptyStart()` above for why such a fixture is worth having and why
 *     it must stay a separate sub-population.
 */
static bool cbcDumpMirFixture(OsiSolverInterface *si, const char *tag = "mir",
  int pass = 0, int options = 0, int inTree = 0, int formulationRows = -1,
  double maxSeconds = 0.0, int level = 0)
{
  if (!si)
    return false;

  const std::string name = cbcMirFixtureName(si);

  if (si->getNumCols() == 0 || si->getNumRows() == 0) {
    printf("[mirfixture] %s.%s: SKIP (empty model: %d rows, %d cols)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols());
    fflush(stdout);
    return false;
  }

  if (!si->isProvenOptimal()) {
    printf("[mirfixture] %s.%s: SKIP (LP not proven optimal; "
           "mixIntRoundPreprocess resolves every range row from getRowActivity() "
           "at CglMixedIntegerRounding2.cpp:508-517 and the resolved sense is what "
           "determineRowType classifies, so without the LP point the row types -- "
           "and hence the whole starting-row list -- are not reproducible)\n",
      name.c_str(), tag);
    fflush(stdout);
    return false;
  }

  CoinWarmStart *rawWs = si->getWarmStart();
  const CoinWarmStartBasis *ws
    = dynamic_cast< const CoinWarmStartBasis * >(rawWs);
  if (!ws) {
    printf("[mirfixture] %s.%s: SKIP (no CoinWarmStartBasis from the solver; MIR2 "
           "does not read basis statuses, but without a basis the replay cannot "
           "reconstruct this LP point and the range-row resolution changes)\n",
      name.c_str(), tag);
    fflush(stdout);
    delete rawWs;
    return false;
  }

  if (ws->getNumStructural() != si->getNumCols()
    || ws->getNumArtificial() != si->getNumRows()) {
    printf("[mirfixture] %s.%s: SKIP (basis dimensions %d/%d do not match the "
           "model's %d cols / %d rows)\n",
      name.c_str(), tag, ws->getNumStructural(), ws->getNumArtificial(),
      si->getNumCols(), si->getNumRows());
    fflush(stdout);
    delete rawWs;
    return false;
  }

  // CBC's configuration, CbcSolverCutSetup.cpp:341-342. CRITERION_ 1 and
  // MULTIPLY_ true -> upperLimit 2; EPSILON_ and UNDEFINED_ are the constructor's
  // (CglMixedIntegerRounding2.cpp:323-324).
  const int criterion = 1;
  const double epsilon = 1.0e-6;
  const int undefined = -1;
  const int upperLimit = 2;

  CbcMirFixtureCounts c;
  cbcMirFixtureCount(si, ws, criterion, epsilon, undefined, upperLimit, c);
  delete rawWs;

  if (c.startingRows == 0 && cbcMirAllowEmptyStart()) {
    // Dumping anyway, as a pure fixed-cost probe. Say so loudly: a downstream reader
    // seeing rowCuts_n == 0 on this fixture must know it is 0 by construction and not
    // because a change suppressed a cut.
    printf("[mirfixture] %s.%s: no starting row, but "
           "CBC_MIR_FIXTURE_ALLOW_EMPTY_START is set, so dumping it as a pure "
           "fixed-cost probe (rowsMix=%d rowsContVB=%d rowsInt=%d: the separation "
           "loop will run zero times and the fixture can produce 0 cuts BY "
           "CONSTRUCTION -- it measures preprocessing and the two O(nnz) matrix "
           "copies only, and cannot support any claim about cut content)\n",
      name.c_str(), tag, c.rowsMix, c.rowsContVB, c.rowsInt);
    fflush(stdout);
  }

  if (c.startingRows == 0 && !cbcMirAllowEmptyStart()) {
    printf("[mirfixture] %s.%s: SKIP (no starting row: %d rows, %d cols, "
           "rowsMix=%d rowsContVB=%d rowsInt=%d, so the outer loop at "
           "CglMixedIntegerRounding2.cpp:866 runs zero times; %d ROW_CONT row(s) "
           "exist but none contains a column with a variable bound, and there are "
           "%d variable upper / %d variable lower bound(s) in total)\n",
      name.c_str(), tag, si->getNumRows(), si->getNumCols(), c.rowsMix,
      c.rowsContVB, c.rowsInt, c.rowsCont, c.vubCount, c.vlbCount);
    fflush(stdout);
    return false;
  }

  double cutoff = COIN_DBL_MAX;
  si->getDblParam(OsiDualObjectiveLimit, cutoff);

  const std::string dir = cbcMirFixtureDir();
  cbcMirFixtureMkdirP(dir);
  const std::string stem = dir + "/" + name + "." + tag;

  // The written file is ".mps.gz": writeMpsNative always gzips and appends the
  // suffix itself.
  int paddedColumns = 0;
  const int mpsRc = cbcMirFixtureWriteMps(si, stem + ".mps", paddedColumns);
  int integerColumns = 0;
  const bool ctOk = cbcMirFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  const int effFormulationRows
    = formulationRows >= 0 ? formulationRows : si->getNumRows();
  const bool metaOk = cbcMirFixtureWriteMeta(si, tag, paddedColumns,
    integerColumns, c, pass, options, inTree, level, effFormulationRows, cutoff,
    maxSeconds, stem + ".meta");
  const bool basOk = cbcMirFixtureWriteBasis(si, stem + ".bas");
  const bool solOk = cbcMirFixtureWriteSol(si, stem + ".sol");

  printf("[mirfixture] %s.%s: %s rows=%d cols=%d int=%d fixed=%d padded=%d "
         "range=%d geFlip=%d mix=%d cont=%d contVB=%d int_=%d varUB=%d varLB=%d "
         "varEQ=%d other=%d undef=%d start=%d vub=%d vlb=%d cmirCalls=%d "
         "knapMax=%d knapSum=%.0f cmirWork=%.0f nnzWork=%.0f rowNzMax=%d "
         "obj=%.15g pass=%d options=%d (mps=%d ctype=%d meta=%d bas=%d sol=%d)\n",
    name.c_str(), tag,
    (mpsRc == 0 && ctOk && metaOk && basOk && solOk) ? "DUMPED" : "PARTIAL",
    si->getNumRows(), si->getNumCols(), integerColumns, c.fixedColumns,
    paddedColumns, c.rangeRows, c.geRowsToFlip, c.rowsMix, c.rowsCont,
    c.rowsContVB, c.rowsInt, c.rowsVarUB, c.rowsVarLB, c.rowsVarEQ, c.rowsOther,
    c.rowsUndefined, c.startingRows, c.vubCount, c.vlbCount, c.cmirCalls,
    c.knapsackIntMax, c.knapsackIntSum, c.cmirWork, c.nnzWork, c.rowNzMax,
    si->getObjValue(), pass, options, mpsRc, (int)ctOk, (int)metaOk, (int)basOk,
    (int)solOk);
  fflush(stdout);

  // The basis is required, not optional -- see the header comment -- so unlike the
  // Probing dump it counts toward success.
  return mpsRc == 0 && ctOk && metaOk && basOk;
}

#endif // CBC_DUMP_MIR_FIXTURE
#endif // CbcMirFixtureDump_H
