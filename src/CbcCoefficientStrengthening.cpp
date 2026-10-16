// Copyright (C) 2026 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcCoefficientStrengthening.hpp"

#include "CoinConflictGraph.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include <unordered_map>
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

// Reference solution loaded by -debugCuts; set before applyLpMethod() so that
// pre-root-LP strengthening can check itself before the OsiRowCutDebugger is
// active. Same globals CbcBoundPropagation.cpp falls back on -- see the
// comment there for why the debugger itself is not yet attached this early
// when driven from the "-debugCuts" CLI switch (it is already active by this
// point when driven through Cbc_activateRowCutDebugger()/mip-debug-cuts).
extern double *debugSolution;
extern int debugNumberColumns;

namespace {

/// Anything at or beyond this magnitude is treated as infinite, matching the
/// convention used throughout Osi/Cgl (see e.g. CglPreProcess).
const double COEFSTR_INFINITY = 1.0e20;

/// Feasibility tolerance. Same role as papilo's Num::feastol: a slack below
/// this is not considered a real slack, and a value this close to its ceiling
/// is snapped up to it.
const double COEFSTR_FEASTOL = 1.0e-6;

/// A column whose bounds are this close together is fixed; tightening its
/// coefficient is exact but pointless, so it is skipped (as papilo does).
const double COEFSTR_EPSILON = 1.0e-9;

/// Skip the clique-cover check (see the class comment) on rows with more
/// unit-coefficient binary columns than this. Measured on 883 corpus
/// instances (mip-sanity-data + a 2017+ MIPLIB/set-packing sample): every
/// instance with a genuine (non-circular) opportunity had its worst row
/// under 20000 such columns; the ones that blew past it (e.g. a single
/// ~194000-column row) drove the greedy cover's wall time into double-digit
/// seconds on their own, all for a single row of one instance.
const size_t CLIQUECOVER_MAX_ROW_VARS = 20000;

/// Deterministic work budget for all clique covers of one run() call, in
/// elements touched (see greedyCliqueCover()). Candidate rows are processed
/// cheapest first, so the budget buys the largest number of rows; once the
/// next row's estimated cost exceeds what is left, the remaining rows keep
/// the plain per-column rule. Counting work instead of wall time makes the
/// strengthened model identical from run to run and machine to machine.
const unsigned long long CLIQUECOVER_WORK_BUDGET = 100000000ULL;

/// One pending right-hand-side change, applied only once the matrix
/// replacement has been verified to have taken effect.
struct RowBoundChange {
  int row;
  double lower;
  double upper;
};

/// Everything the per-row analysis produces; shared by the clique-cover
/// planning pass and the strengthening pass so both see identical numbers.
struct RowInfo {
  bool hasUpper;
  double scale;
  double rhs;
  double slack;
  double minActivityOthers;
  bool minOthersUnbounded;
  std::vector< int > unitBinaryCols;
};

/*! \brief Greedy clique COVER of \p remaining (conflict-graph node ids,
 *         ascending, unique -- ownership taken and consumed).
 *
 * Repeatedly grows one maximal clique (seed = the smallest remaining node,
 * then keep intersecting the candidate set with each newly-added member's
 * neighbourhood, per CoinConflictGraph::conflictingNodes()), removes it from
 * \p remaining, and repeats. Each step only shrinks the sorted candidate set
 * via std::set_intersection, so the total work is bounded by the sum of
 * neighbourhood sizes actually touched rather than the O(k^2) of checking
 * every pair directly.
 *
 * Both stopping rules are deterministic (they count work, never time):
 *  - \p workLimit: once the work counter reaches it, every still-uncovered
 *    node is conservatively counted as its own singleton clique (always a
 *    valid, if not maximally tight, upper bound on how many of them can be
 *    simultaneously 1).
 *  - \p maxUseful: the largest cover size that would still change anything
 *    for the caller. A cover of more than this many cliques is useless, and
 *    the count only grows, so the moment it is certain to exceed it the
 *    search stops with \p pruned set.
 *
 * \param work  Incremented by the elements touched (neighbour lists
 *              retrieved, candidate sets intersected, remaining set copied).
 * \return Number of cliques used to cover \p remaining (meaningless if
 *         \p pruned).
 */
size_t greedyCliqueCover(const CoinConflictGraph *cg,
  std::vector< size_t > remaining, size_t *tempBuf, char *ivBuf,
  unsigned long long workLimit, unsigned long long &work, size_t maxUseful,
  bool &pruned)
{
  size_t nCliques = 0;
  pruned = false;
  while (!remaining.empty()) {
    if (nCliques + 1 > maxUseful) {
      pruned = true;
      return nCliques + 1;
    }
    if (work >= workLimit) {
      nCliques += remaining.size();
      break;
    }

    const size_t seed = remaining.front();
    std::pair< size_t, const size_t * > seedNb = cg->conflictingNodes(seed, tempBuf, ivBuf);
    work += seedNb.first + remaining.size();

    std::vector< size_t > candidates;
    candidates.reserve(remaining.size());
    std::set_intersection(remaining.begin() + 1, remaining.end(),
      seedNb.second, seedNb.second + seedNb.first,
      std::back_inserter(candidates));

    std::vector< size_t > cliqueMembers;
    cliqueMembers.push_back(seed);
    while (!candidates.empty() && work < workLimit) {
      const size_t next = candidates.front();
      cliqueMembers.push_back(next);

      std::pair< size_t, const size_t * > nextNb = cg->conflictingNodes(next, tempBuf, ivBuf);
      work += nextNb.first + candidates.size();
      std::vector< size_t > narrowed;
      narrowed.reserve(candidates.size());
      std::set_intersection(candidates.begin() + 1, candidates.end(),
        nextNb.second, nextNb.second + nextNb.first,
        std::back_inserter(narrowed));
      candidates.swap(narrowed);
    }
    nCliques++;

    // cliqueMembers is built in ascending order (each pick is the smallest
    // surviving candidate of a shrinking sorted set), so a straight
    // set_difference removes it from `remaining` without re-sorting.
    std::vector< size_t > nextRemaining;
    nextRemaining.reserve(remaining.size() - cliqueMembers.size());
    std::set_difference(remaining.begin(), remaining.end(),
      cliqueMembers.begin(), cliqueMembers.end(),
      std::back_inserter(nextRemaining));
    remaining.swap(nextRemaining);
  }
  return nCliques;
}

} // end anonymous namespace

CbcCoefficientStrengthening::CbcCoefficientStrengthening()
  : nCoefficients_(0)
  , nRows_(0)
  , timeUsed_(0.0)
  , nCliqueCoverRows_(0)
  , cliqueCoverReduction_(0)
  , cliqueCoverTimeUsed_(0.0)
  , cliqueCoverCandidates_(0)
  , cliqueCoverPruned_(0)
  , cliqueCoverBudgetSkipped_(0)
  , cliqueCoverWork_(0)
{
}

bool CbcCoefficientStrengthening::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel)
{
  nCoefficients_ = 0;
  nRows_ = 0;
  timeUsed_ = 0.0;
  nCliqueCoverRows_ = 0;
  cliqueCoverReduction_ = 0;
  cliqueCoverTimeUsed_ = 0.0;
  cliqueCoverCandidates_ = 0;
  cliqueCoverPruned_ = 0;
  cliqueCoverBudgetSkipped_ = 0;
  cliqueCoverWork_ = 0;

  if (!solver || solver->getNumRows() == 0 || solver->getNumCols() == 0)
    return false;

  const double t0 = CoinWallclockTime();

  const int numberRows = solver->getNumRows();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();

  // Scan the solver's own row copy read-only, and only take a mutable copy of
  // the matrix once a change is actually known to be needed. On this corpus
  // the pass finds nothing on 387 of 471 instances, and copying the matrix
  // eagerly made those instances 57% of its total cost -- dt_optimization
  // (1.0M nonzeros, no change) spent 5.7 ms of its 5.9 ms purely on a copy
  // that was then discarded. The scan itself is cheap: it reads bounds and
  // coefficients and writes nothing.
  const CoinPackedMatrix *rowCopy = solver->getMatrixByRow();
  const double *element = rowCopy->getElements();
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();

  // Deferred mutable copy plus the writes destined for it. A change is
  // recorded here during the scan and replayed into the copy afterwards,
  // which keeps the copy off the no-change path entirely.
  //
  // Note this must be a copy of the *pre-change* matrix, so it has to be
  // taken before any element is written -- hence recording the writes rather
  // than applying them as they are found.
  struct ElementChange {
    CoinBigIndex position;
    double value;
  };
  std::vector< ElementChange > elementChanges;

  std::vector< RowBoundChange > rowBoundChanges;
  int firstRow = -1, firstColumn = -1;
  double firstValue = 0.0;

  // Conflict graph for the clique-cover check (see the class comment). Only
  // ever *read* here -- built earlier in the pre-root-LP phase by clique
  // merging (when -clqStrengthening is before/both) -- and null otherwise, in
  // which case every row below just falls back to the plain per-column rule.
  const CoinConflictGraph *cgraph = solver->getCGraph();
  std::vector< size_t > cliqueTempBuf;
  std::vector< char > cliqueIvBuf;
  if (cgraph) {
    cliqueTempBuf.resize(cgraph->size());
    cliqueIvBuf.assign(cgraph->size(), 0);
  }

  auto analyzeRow = [&](int iRow, RowInfo &ri) -> bool {
    const CoinBigIndex start = rowStart[iRow];
    const int length = rowLength[iRow];
    // A singleton row carries no slack argument: shrinking its only
    // coefficient would just be a bound change, which is bound propagation's
    // job (and it has already run).
    if (length <= 1)
      return false;

    // Exactly one finite side. A range row has two slacks and neither of them
    // alone justifies the change; a free row has none.
    const bool hasUpper = rowUpper[iRow] < COEFSTR_INFINITY;
    const bool hasLower = rowLower[iRow] > -COEFSTR_INFINITY;
    if (hasUpper == hasLower)
      return false;

    // Normalise to a.x <= rhs: a >= row is negated wholesale.
    const double scale = hasUpper ? 1.0 : -1.0;
    double rhs = hasUpper ? rowUpper[iRow] : -rowLower[iRow];

    // Maximum activity. An infinite contribution makes the row unbounded
    // above, so there is no slack to work with at all. Also collects:
    //  - unitBinaryCols: binary columns with coefficient exactly +1 (after
    //    normalisation) -- candidates for the clique-cover check below;
    //  - minActivityOthers / minOthersUnbounded: the row's *minimum* possible
    //    activity contributed by every OTHER (non-unitBinaryCols) column,
    //    i.e. each held at whichever bound frees up the most room for the
    //    row to stay feasible. This is the anti-circularity guard's input
    //    (see the comment at its use below) and is unrelated to maxActivity,
    //    which uses the opposite (most-favorable-for-violation) bound.
    double maxActivity = 0.0;
    double minActivityOthers = 0.0;
    bool infiniteActivity = false;
    bool minOthersUnbounded = false;
    std::vector< int > unitBinaryCols;
    for (CoinBigIndex j = start; j < start + length; j++) {
      const double value = element[j] * scale;
      const int iColumn = column[j];
      if (value > 0.0) {
        if (colUpper[iColumn] >= COEFSTR_INFINITY) {
          infiniteActivity = true;
          break;
        }
        maxActivity += value * colUpper[iColumn];
        // A column already fixed at 0 contributes nothing to maxActivity
        // (colUpper[iColumn] == 0 above), yet the clique-cover reduction
        // below subtracts (k - p) assuming every one of the k members
        // contributes exactly 1 to the naive maxActivity. Counting a
        // fixed-at-0 member in k would inflate that reduction with no
        // matching maxActivity to justify it -- so such columns are excluded
        // from the candidate set entirely (they can never be part of the
        // "at most p of these can be 1" argument's k side, since they can
        // never be 1 at all).
        const bool isUnitBinary = cgraph != NULL
          && std::fabs(value - 1.0) <= COEFSTR_FEASTOL
          && solver->isBinary(iColumn)
          && colUpper[iColumn] >= 1.0 - COEFSTR_FEASTOL;
        if (isUnitBinary) {
          unitBinaryCols.push_back(iColumn);
        } else if (colLower[iColumn] <= -COEFSTR_INFINITY) {
          minOthersUnbounded = true;
        } else {
          minActivityOthers += value * colLower[iColumn];
        }
      } else if (value < 0.0) {
        if (colLower[iColumn] <= -COEFSTR_INFINITY) {
          infiniteActivity = true;
          break;
        }
        maxActivity += value * colLower[iColumn];
        if (colUpper[iColumn] >= COEFSTR_INFINITY)
          minOthersUnbounded = true;
        else
          minActivityOthers += value * colUpper[iColumn];
      }
    }
    if (infiniteActivity)
      return false;

    // The slack is what every coefficient in this row is capped at. It is
    // invariant under the change itself (the coefficient loses exactly what
    // the rhs loses), so it is computed once and the order of the changes
    // below does not matter.
    double slack = maxActivity - rhs;
    if (slack <= COEFSTR_FEASTOL)
      return false;
    ri.hasUpper = hasUpper;
    ri.scale = scale;
    ri.rhs = rhs;
    ri.slack = slack;
    ri.minActivityOthers = minActivityOthers;
    ri.minOthersUnbounded = minOthersUnbounded;
    ri.unitBinaryCols.swap(unitBinaryCols);
    return true;
  };

  // Clique-cover planning (see the class comment's "Conflict-aware slack"
  // section). The cover only ever lowers the slack used for the row's OTHER
  // integer columns, so most rows cannot benefit at all; for those that can,
  // the cover is worth running only if it ends with few enough cliques, which
  // gives an early-abort threshold. Candidates are then covered cheapest
  // first under a deterministic work budget.
  std::unordered_map< int, size_t > coverSize;
  if (cgraph) {
    struct Candidate {
      int row;
      unsigned long long cost;
      size_t maxUseful;
    };
    std::vector< Candidate > candidates;
    RowInfo ri;
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (!analyzeRow(iRow, ri))
        continue;
      const size_t k = ri.unitBinaryCols.size();
      if (k < 2 || k > CLIQUECOVER_MAX_ROW_VARS)
        continue;
      const double cap = ri.minOthersUnbounded
        ? COEFSTR_INFINITY
        : (ri.rhs - ri.minActivityOthers);
      if (cap < 2.0 - COEFSTR_FEASTOL)
        continue;

      std::sort(ri.unitBinaryCols.begin(), ri.unitBinaryCols.end());
      double largestOther = 0.0;
      const CoinBigIndex start = rowStart[iRow];
      for (CoinBigIndex j = start; j < start + rowLength[iRow]; j++) {
        const int iColumn = column[j];
        if (!solver->isInteger(iColumn)
          || colLower[iColumn] <= -COEFSTR_INFINITY
          || colUpper[iColumn] >= COEFSTR_INFINITY
          || std::binary_search(ri.unitBinaryCols.begin(),
               ri.unitBinaryCols.end(), iColumn))
          continue;
        if (std::floor(colUpper[iColumn] + COEFSTR_FEASTOL)
            - std::ceil(colLower[iColumn] - COEFSTR_FEASTOL)
          < 1.0 - COEFSTR_FEASTOL)
          continue;
        largestOther = std::max(largestOther, std::fabs(element[j]));
      }
      // A cover of p cliques gives slack' = slack - (k - p); it changes
      // something only if some other coefficient exceeds slack', i.e.
      // p < largestOther - slack + k. The row's own guard needs
      // (k - p) < slack, i.e. p > k - slack.
      const double upperP = largestOther - ri.slack + static_cast< double >(k);
      const double lowerP = static_cast< double >(k) - ri.slack;
      if (upperP <= 1.0 + COEFSTR_FEASTOL)
        continue;
      size_t maxUseful = static_cast< size_t >(std::ceil(upperP - COEFSTR_FEASTOL)) - 1;
      maxUseful = std::min(maxUseful, k - 1);
      if (static_cast< double >(maxUseful) <= lowerP + COEFSTR_FEASTOL)
        continue;

      Candidate c;
      c.row = iRow;
      c.maxUseful = maxUseful;
      c.cost = 0;
      for (size_t i = 0; i < k; i++)
        c.cost += cgraph->degree(static_cast< size_t >(ri.unitBinaryCols[i]));
      candidates.push_back(c);
    }
    cliqueCoverCandidates_ = static_cast< int >(candidates.size());

    std::sort(candidates.begin(), candidates.end(),
      [](const Candidate &a, const Candidate &b) {
        return a.cost != b.cost ? a.cost < b.cost : a.row < b.row;
      });

    const unsigned long long budget = CLIQUECOVER_WORK_BUDGET;
    for (size_t ic = 0; ic < candidates.size(); ic++) {
      const Candidate &c = candidates[ic];
      if (cliqueCoverWork_ >= budget) {
        cliqueCoverBudgetSkipped_ = static_cast< int >(candidates.size() - ic);
        break;
      }
      if (!analyzeRow(c.row, ri))
        continue;
      std::vector< size_t > nodes(ri.unitBinaryCols.begin(), ri.unitBinaryCols.end());
      std::sort(nodes.begin(), nodes.end());
      bool pruned = false;
      const double coverT0 = CoinWallclockTime();
      const size_t p = greedyCliqueCover(cgraph, nodes, cliqueTempBuf.data(),
        cliqueIvBuf.data(), budget, cliqueCoverWork_, c.maxUseful, pruned);
      cliqueCoverTimeUsed_ += CoinWallclockTime() - coverT0;
      if (pruned)
        cliqueCoverPruned_++;
      else
        coverSize[c.row] = p;
    }
  }

  for (int iRow = 0; iRow < numberRows; iRow++) {
    const CoinBigIndex start = rowStart[iRow];
    const int length = rowLength[iRow];
    RowInfo rowInfo;
    if (!analyzeRow(iRow, rowInfo))
      continue;
    const bool hasUpper = rowInfo.hasUpper;
    const double scale = rowInfo.scale;
    double rhs = rowInfo.rhs;
    double slack = rowInfo.slack;
    std::vector< int > &unitBinaryCols = rowInfo.unitBinaryCols;

    // Slack used for unitBinaryCols themselves -- always the naive
    // (non-clique-adjusted) value; see the comment where it is consumed
    // below for why.
    double slackForUnitBinaryCols = slack;

    // Clique-cover result for this row (computed by the planning pass above).
    const std::unordered_map< int, size_t >::const_iterator coverIt
      = coverSize.find(iRow);
    if (coverIt != coverSize.end()) {
      const int k = static_cast< int >(unitBinaryCols.size());
      const int p = static_cast< int >(coverIt->second);
      // p < k means at most p (not k) of these columns can be 1 at once,
      // so the row's true maximum activity is maxActivity - (k - p), i.e.
      // slack shrinks by the same amount. Guard against the (never
      // observed on the corpus, but not provably impossible) case where
      // that would drive slack to <= 0: unlike the naive maxActivity bound,
      // a clique cover is derived from facts *external* to this row, so a
      // non-positive result here would mean the row is actually redundant
      // once those external constraints are accounted for -- a real but
      // different finding (row redundancy) that the per-column rule below
      // is not equipped to act on safely, so it is left for another pass.
      //
      // This adjusted `slack` must ONLY be used to shrink OTHER (non-
      // unitBinaryCols) coefficients in the row below, e.g. a companion
      // big-M column -- never the unitBinaryCols themselves. The bound
      // "at most p of these k can be 1 simultaneously" justifies removing
      // (k - p) from the row's rhs *once*, in aggregate; it does not
      // justify reducing every one of the k members' own coefficients
      // independently by the same amount, since the standard per-column
      // shrink rule below decrements the rhs again for *each* column it
      // touches (valid when columns are independent, but these k are
      // mutually exclusive, not independent -- applying it to all of them
      // compounds the rhs reduction k-fold instead of the single time the
      // clique fact actually buys, incorrectly cutting off feasible
      // points, e.g. an all-clique-members-zero point at a high value of
      // some other row column). See CBC_DISABLE_CLIQUECOVER's history /
      // the physiciansched3-3 investigation for a worked counterexample.
      if (p < k && (k - p) < slack - COEFSTR_FEASTOL) {
        slack -= (k - p);
        nCliqueCoverRows_++;
        cliqueCoverReduction_ += (k - p);
      }
    }

    // Snap up to an integral slack when we are within tolerance of one, so a
    // slack that is really 3 but computed as 2.9999995 does not shave a
    // sliver off the feasible set. Rounding is only ever *upwards*: a larger
    // cap is a weaker (safe) tightening, a smaller one can cut off integer
    // points. Non-integral slacks are left alone -- the rule does not need
    // them to be integral.
    const double slackCeil = std::ceil(slack);
    if (std::fabs(slackCeil - slack) <= COEFSTR_FEASTOL)
      slack = slackCeil;
    const double slackForUnitBinaryColsCeil = std::ceil(slackForUnitBinaryCols);
    if (std::fabs(slackForUnitBinaryColsCeil - slackForUnitBinaryCols) <= COEFSTR_FEASTOL)
      slackForUnitBinaryCols = slackForUnitBinaryColsCeil;

    std::vector< int > unitBinaryColsSorted(unitBinaryCols.begin(), unitBinaryCols.end());
    std::sort(unitBinaryColsSorted.begin(), unitBinaryColsSorted.end());

    bool changedRow = false;
    for (CoinBigIndex j = start; j < start + length; j++) {
      const int iColumn = column[j];
      if (!solver->isInteger(iColumn))
        continue;
      // Both bounds must be finite: the rhs adjustment is stated in terms of
      // the bound the column attains at maximum activity.
      if (colLower[iColumn] <= -COEFSTR_INFINITY || colUpper[iColumn] >= COEFSTR_INFINITY)
        continue;

      // The adjustment is stated at an *attainable integer* bound: the
      // argument is that one step in from it the row is already redundant,
      // and a step is only a step if the bound is integral. Bound
      // propagation normally leaves integral bounds on integer columns, but a
      // bound of 0.5 would otherwise give a step of 0.5 and shave a sliver
      // off the feasible set, so round inwards and require a real step to
      // exist -- which also excludes fixed columns, where shrinking the
      // coefficient buys nothing.
      const double lower = std::ceil(colLower[iColumn] - COEFSTR_FEASTOL);
      const double upper = std::floor(colUpper[iColumn] + COEFSTR_FEASTOL);
      if (upper - lower < 1.0 - COEFSTR_FEASTOL)
        continue;

      const double value = element[j] * scale;
      // unitBinaryCols members are mutually exclusive (per the clique-cover
      // fact used above, when it applied), not independent, so they must use
      // the un-adjusted slack -- see the long comment above.
      const bool isUnitBinaryCol = std::binary_search(
        unitBinaryColsSorted.begin(), unitBinaryColsSorted.end(), iColumn);
      const double effectiveSlack = isUnitBinaryCol ? slackForUnitBinaryCols : slack;
      if (std::fabs(value) <= effectiveSlack + COEFSTR_FEASTOL)
        continue;

      const double newValue = (value > 0.0 ? effectiveSlack : -effectiveSlack);
      // At the bound the column takes at maximum activity the two rows agree
      // exactly; one integer step away from it both are redundant.
      rhs -= (value - newValue) * (value > 0.0 ? upper : lower);
      // Recorded, not written: `element` still points at the solver's own row
      // copy. Safe to defer because the slack is invariant under the change
      // and each coefficient is considered exactly once, so nothing later in
      // the scan depends on this write having landed.
      ElementChange elementChange;
      elementChange.position = j;
      elementChange.value = newValue * scale;
      elementChanges.push_back(elementChange);

      if (logLevel >= 3)
        printf("  Coefficient strengthening: row %d col %d: %g -> %g "
               "(slack %g, rhs now %g)\n",
          iRow, iColumn, value, newValue, effectiveSlack, rhs);

      if (firstRow < 0) {
        firstRow = iRow;
        firstColumn = iColumn;
        firstValue = elementChange.value;
      }
      nCoefficients_++;
      changedRow = true;
    }

    if (changedRow) {
      nRows_++;
      RowBoundChange change;
      change.row = iRow;
      change.lower = hasUpper ? rowLower[iRow] : -rhs;
      change.upper = hasUpper ? rhs : rowUpper[iRow];
      rowBoundChanges.push_back(change);
    }
  }

  if (!nCoefficients_) {
    timeUsed_ = CoinWallclockTime() - t0;
    return false;
  }

  // Only now is the copy worth its cost. Poking the single entries with
  // modifyCoefficient() instead would leave the solver's cached row copy stale
  // -- ClpModel::modifyCoefficient() only clears whatsChanged_ and never calls
  // freeCachedResults() -- so a coefficient written that way reads back with
  // its old value through getMatrixByRow(). replaceMatrix() does free the
  // caches, which is also why Cgl's own tightener batches this way.
  CoinPackedMatrix copy = *rowCopy;
  double *mutableElement = copy.getMutableElements();
  for (size_t i = 0; i < elementChanges.size(); i++)
    mutableElement[elementChanges[i].position] = elementChanges[i].value;

  solver->replaceMatrixOptional(copy);

  // replaceMatrixOptional() is a no-op in the OsiSolverInterface base class,
  // so for a solver that does not implement it the matrix would silently keep
  // its original coefficients. Adjusting the row bounds anyway would then cut
  // off feasible points, so confirm the replacement landed before touching
  // them, and leave the model exactly as it was if it did not.
  bool applied = false;
  {
    const CoinPackedMatrix *check = solver->getMatrixByRow();
    if (check && firstRow >= 0 && firstRow < check->getNumRows()) {
      const CoinBigIndex start = check->getVectorStarts()[firstRow];
      const int length = check->getVectorLengths()[firstRow];
      const int *indices = check->getIndices();
      const double *values = check->getElements();
      for (CoinBigIndex j = start; j < start + length; j++) {
        if (indices[j] == firstColumn) {
          applied = (std::fabs(values[j] - firstValue) < COEFSTR_EPSILON);
          break;
        }
      }
    }
  }

  if (!applied) {
    if (logLevel >= 1)
      printf("  Coefficient strengthening: solver did not accept the modified "
             "matrix — %d changes discarded.\n",
        nCoefficients_);
    nCoefficients_ = 0;
    nRows_ = 0;
    timeUsed_ = CoinWallclockTime() - t0;
    return false;
  }

  for (size_t i = 0; i < rowBoundChanges.size(); i++) {
    const RowBoundChange &change = rowBoundChanges[i];
    if (change.lower > -COEFSTR_INFINITY)
      solver->setRowLower(change.row, change.lower);
    if (change.upper < COEFSTR_INFINITY)
      solver->setRowUpper(change.row, change.upper);
  }

  // Debugger validation: every changed row must still contain the reference
  // solution once evaluated with its NEW coefficients against its NEW
  // bounds. This is a purely algebraic invariant of the rule above (see the
  // class comment) -- unlike a cut, nothing here depends on which subtree of
  // the B&B tree we are in, so it is checked unconditionally whenever a
  // debugger/-debugCuts solution is available, not just on the optimal path.
  //
  // Fall back to the debugSolution global (see the comment above the
  // extern declarations) when no OsiRowCutDebugger has been attached to
  // `solver` yet -- this pass runs ahead of the root LP solve, so a
  // "-debugCuts" driven debugger is not attached this early, while
  // Cbc_activateRowCutDebugger()/mip-debug-cuts attach it before Cbc_solve()
  // is even called, so a debugger already exists there.
  //
  // Deliberately use getRowCutDebuggerAlways() and NOT getRowCutDebugger():
  // the latter gates on OsiRowCutDebugger::onOptimalPath(), which only
  // checks INTEGER columns with a loose 1e-3 tolerance and has been observed
  // to be unreliable (can return false even though the current bounds are
  // still fully consistent with the reference solution) -- see the matching
  // fixes/comments in CbcModel.cpp (reducedCostFix(), fixFromGlobalCuts())
  // and CbcBoundPropagation::run(). Using the unreliable check here would
  // silently skip this validation whenever it spuriously returns false,
  // which is exactly the failure mode this diagnostic exists to catch.
  const OsiRowCutDebugger *debugger = solver->getRowCutDebuggerAlways();
  const double *optSol = debugger
    ? debugger->optimalSolution()
    : (debugSolution && debugNumberColumns == solver->getNumCols()
          ? debugSolution
          : nullptr);
  if (optSol) {
    const double *checkElement = copy.getElements();
    const int *checkColumn = copy.getIndices();
    const CoinBigIndex *checkRowStart = copy.getVectorStarts();
    const int *checkRowLength = copy.getVectorLengths();
    const double tol = 1.0e-6;
    for (size_t i = 0; i < rowBoundChanges.size(); i++) {
      const RowBoundChange &change = rowBoundChanges[i];
      const int iRow = change.row;
      double activity = 0.0;
      const CoinBigIndex start = checkRowStart[iRow];
      const int length = checkRowLength[iRow];
      for (CoinBigIndex j = start; j < start + length; j++)
        activity += checkElement[j] * optSol[checkColumn[j]];
      const bool violatesUpper = change.upper < COEFSTR_INFINITY && activity > change.upper + tol;
      const bool violatesLower = change.lower > -COEFSTR_INFINITY && activity < change.lower - tol;
      if (violatesUpper || violatesLower) {
        printf("coefStrengthening BAD ROW: row %d activity=%.12g at reference "
               "solution but new bounds are [%.12g,%.12g] (violates %s)\n",
          iRow, activity, change.lower, change.upper,
          violatesUpper ? "upper" : "lower");
        fflush(stdout);
      }
    }
  }

  timeUsed_ = CoinWallclockTime() - t0;
  return true;
}
