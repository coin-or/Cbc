/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * @file CbcCutPoolFilter.cpp
 * @brief Generic post-generation cut-pool filter for CbcCutGenerator
 *
 **/

#include "CbcCutPoolFilter.hpp"

#include "CoinCutPool.hpp"
#include "CoinPackedVector.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"

#include <cstdlib>
#include <vector>

namespace {

const double CBC_CUTPOOL_FILTER_INF = 1.0e30;

int envInt(const char *name, int def)
{
  const char *v = getenv(name);
  return v ? atoi(v) : def;
}

double envDouble(const char *name, double def)
{
  const char *v = getenv(name);
  return v ? atof(v) : def;
}

} // namespace

bool cbcFilterGeneratedCuts(OsiCuts &cs, int firstRowCut, const double *x,
  int numCols, const char *generatorTag)
{
  const int lastRowCut = cs.sizeRowCuts();
  const int nCandidates = lastRowCut - firstRowCut;
  if (nCandidates <= 0)
    return false;

  // Same env-var naming pattern (and rationale) as CglBKClique's
  // CBC_CLIQUE_POOL_* precedent: small models / small candidate counts are
  // exempt, since filtering only pays for itself when there are many
  // candidates competing for a handful of "best cut for this column" slots.
  static const int minCols = envInt("CBC_CUTPOOL_FILTER_MIN_COLS", 500);
  static const int minCandidates = envInt("CBC_CUTPOOL_FILTER_MIN_CANDIDATES", 20);
  static const bool alwaysFilter = envInt("CBC_CUTPOOL_FILTER_ALWAYS", 0) != 0;
  // Unlike CglBKClique's own clique-cut parallelism filter (disabled by
  // default -- a 442-instance sweep found no net win for clique cuts
  // specifically), the 2026-09 mip-sanity-data sweep for these four
  // generators found MAX_PARALLELISM=0.7 gave the single best combined
  // primal-gap/efficiency result of every variant tried (see
  // ROOT-FIXTURES.md's "Cut-pool filtering" section), so it ships as the
  // default here rather than left disabled pending further confirmation.
  static const double maxParallelism = envDouble("CBC_CUTPOOL_FILTER_MAX_PARALLELISM", 0.7);

  const bool smallModel = numCols < minCols;
  if (smallModel && !alwaysFilter)
    return false;

  // First pass: classify which candidates fit the pool's single-sided
  // "ax <= rhs" model. A cut with only ub() binding (lb() effectively -inf)
  // is used as-is; a cut with only lb() binding (ub() effectively +inf) is
  // negated to "-ax <= -lb" form for scoring purposes. Genuine ranged or
  // equality cuts (both bounds finite and used) don't fit this model and
  // are left untouched -- unconditionally kept.
  std::vector<signed char> sign(nCandidates, 0);
  int nFilterable = 0;
  for (int k = firstRowCut; k < lastRowCut; k++) {
    const OsiRowCut *rc = cs.rowCutPtr(k);
    const double lb = rc->lb();
    const double ub = rc->ub();
    const int idx = k - firstRowCut;
    if (rc->row().getNumElements() == 0)
      continue;
    if (ub < CBC_CUTPOOL_FILTER_INF && lb <= -CBC_CUTPOOL_FILTER_INF) {
      sign[idx] = 1;
      nFilterable++;
    } else if (lb > -CBC_CUTPOOL_FILTER_INF && ub >= CBC_CUTPOOL_FILTER_INF) {
      sign[idx] = -1;
      nFilterable++;
    }
  }

  if (nFilterable == 0)
    return false;
  if (nFilterable < minCandidates && !alwaysFilter)
    return false;

  CoinCutPool pool(x, numCols, generatorTag);
  pool.setMaxParallelism(maxParallelism);

  std::vector<int> idxBuf(numCols);
  std::vector<double> coefBuf(numCols);
  for (int k = firstRowCut; k < lastRowCut; k++) {
    const int idx = k - firstRowCut;
    if (!sign[idx])
      continue;
    const OsiRowCut *rc = cs.rowCutPtr(k);
    const CoinPackedVector &row = rc->row();
    const int nz = row.getNumElements();
    const int *idxs = row.getIndices();
    const double *els = row.getElements();
    const double s = sign[idx];
    const double rhs = (s > 0.0) ? rc->ub() : -rc->lb();
    for (int t = 0; t < nz; t++) {
      idxBuf[t] = idxs[t];
      coefBuf[t] = s * els[t];
    }
    pool.add(idxBuf.data(), coefBuf.data(), nz, rhs, static_cast< size_t >(idx));
  }

  pool.removeNullCuts();
  pool.filterByParallelism();

  std::vector< bool > keep(nCandidates, false);
  for (size_t i = 0; i < pool.numCuts(); i++)
    keep[pool.cutTag(i)] = true;

  bool removedAny = false;
  // Erase from highest to lowest index so earlier indices stay valid
  // (OsiCuts::eraseRowCut() shifts subsequent entries down).
  for (int k = lastRowCut - 1; k >= firstRowCut; k--) {
    const int idx = k - firstRowCut;
    if (sign[idx] && !keep[idx]) {
      cs.eraseRowCut(k);
      removedAny = true;
    }
  }
  return removedAny;
}
