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
  int numCols, int numElements, const char *generatorTag)
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
  // Secondary "cheap to reoptimize" signal: constraint-matrix nonzero
  // count. Derived offline (2026-09) by logging every root cut-generation
  // round's actual LP reoptimization CPU time
  // (CbcModel::lastCutRoundResolveTime(), surfaced via
  // CBC_LOG_ROOT_RESOLVE_TIME) across the 442-instance mip-sanity-data
  // set and comparing against numElements: instances with
  // numElements<=8000 show a p90 reoptimization cost (0.032s) matching
  // the p90 already accepted for cols<minCols instances (0.034s), and a
  // worst case (0.16s) well under that group's worst case (0.72s) -- so
  // this class of instances is genuinely as cheap to reoptimize as the
  // ones already exempted by the cols<500 gate. However, a controlled
  // apples-to-apples validation sweep of *actually skipping filtering*
  // there (not just the reoptimization-time proxy) found the effect is a
  // wash once B&B-path noise is averaged out: a single-seed run showed a
  // spread-looking "32 improved / 22 worse" split on this zone's 93
  // instances, but re-running with 5 independent -randomSeed repeats per
  // config and averaging per-instance (cutfilter-sweep --repeats=5; see
  // ROOT-FIXTURES.md) collapsed that entirely -- 0.00pp net dual/primal
  // gap-closed delta, and a negligible +0.05s mean bbTime change -- i.e.
  // the earlier mixed signal was branching tie-break noise from cut
  // content changing node visitation order, not a real effect either
  // way. Ships OFF (0, i.e. no secondary gate) since no benefit was
  // demonstrated (the harmlessness is not in question, but there's no
  // reason to add the extra gate/complexity for zero measured gain); the
  // env var lets this be revisited if a different instance mix or a
  // harder validation set ever shows a real signal. This is intentionally
  // a purely static/offline-derived threshold -- CPU timing is never read
  // back within a live solve to make this decision (see
  // lastCutRoundResolveTime()'s doc comment for why that would be
  // unsound: reoptimization CPU time is noisy and machine-load-dependent,
  // so gating on it live would make cut selection non-reproducible).
  static const int minElements = envInt("CBC_CUTPOOL_FILTER_MIN_ELEMENTS", 0);
  // 2026-09/10 extensive parameter sweep (cutfilter-extensive-r3, 19
  // variants x 3 -randomSeed repeats; then re-confirmed with 5 repeats in
  // cutfilter-combo-r5) found MIN_CANDIDATES=10 (down from the initial
  // CglBKClique-parity default of 20) dominates: dual gap closed +0.30pp,
  // primal gap closed +0.33pp, mean bbTime *improved* (-0.04s) -- a
  // strictly better, essentially free change (relaxing this gate lets
  // filtering kick in for more, somewhat smaller cut-candidate rounds,
  // which turns out to pay for itself in this instance mix). Ships as the
  // new default.
  static const int minCandidates = envInt("CBC_CUTPOOL_FILTER_MIN_CANDIDATES", 10);
  static const bool alwaysFilter = envInt("CBC_CUTPOOL_FILTER_ALWAYS", 0) != 0;
  // Same 2026-09/10 sweep found MAX_PARALLELISM has a real, repeat-averaged
  // monotonic trend across 0.3-1.0: higher (less aggressive orthogonality
  // filtering) consistently improves dual gap closed at a small bbTime
  // cost. 0.9 gave the best combined result together with
  // MIN_CANDIDATES=10 above (cand10_par09 in cutfilter-combo-r5, 5
  // repeats): dual +0.22pp, primal +0.70pp, bbTime +0.03s -- dominates the
  // previous 0.7 default (itself found in an earlier, less thorough
  // sweep) on both bound metrics for a negligible time cost. Ships as the
  // new default, superseding the previous 0.7.
  static const double maxParallelism = envDouble("CBC_CUTPOOL_FILTER_MAX_PARALLELISM", 0.9);

  const bool smallModel = numCols < minCols || numElements <= minElements;
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
