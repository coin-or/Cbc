/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Generalizes the CoinCutPool-based fitness filtering that CglBKClique
 * already applies to clique cuts (see Cgl/src/CglBKClique/CglBKClique.cpp,
 * insertCuts()) to any row-cut generator that can flood a single round with
 * many cuts of similar/redundant support -- currently used for CglGomory,
 * CglMixedIntegerRounding2, CglTwomir and CglProbing (see
 * CbcCutGenerator::generateCuts()).
 *
 * Unlike CglBKClique (which builds CoinCuts directly and re-emits fresh
 * OsiRowCuts from the pool), this filter operates on cuts that already
 * exist as OsiRowCut objects appended to an OsiCuts collection: it feeds
 * each candidate's row into a CoinCutPool (tagging it with its original
 * position), lets the pool's best-fitness-per-column contest (and,
 * optionally, its parallelism filter) decide which candidates to keep,
 * then erases the losers directly from the OsiCuts collection -- so
 * surviving cuts keep their original OsiRowCut identity/metadata
 * (effectiveness, globally-valid flag, etc.) untouched.
 *
 * @file CbcCutPoolFilter.hpp
 * @brief Generic post-generation cut-pool filter for CbcCutGenerator
 *
 **/

#ifndef CBCCUTPOOLFILTER_HPP
#define CBCCUTPOOLFILTER_HPP

class OsiCuts;

/**
 * Filter row cuts appended to `cs` at index >= firstRowCut (i.e. the ones
 * just produced by one generateCuts() call) using a CoinCutPool, keeping
 * only cuts that win the pool's best-fitness-per-column contest (and,
 * optionally, its parallelism/orthogonality filter). Cuts that don't fit
 * the pool's single-sided "ax <= rhs, violated at x" model (genuine ranged
 * or equality cuts) are left untouched -- only >=/<= cuts are candidates
 * for removal.
 *
 * Gated off (a no-op) for small models and small candidate counts, mirroring
 * CglBKClique's precedent -- filtering only pays for itself when there are
 * many candidates to choose among. See CbcCutPoolFilter.cpp for the exact
 * thresholds and their CBC_CUTPOOL_FILTER_* env var overrides.
 *
 * @param cs cuts collection to filter in place.
 * @param firstRowCut first index (inclusive) of the row cuts to consider;
 *        cuts before this index are left untouched.
 * @param x current LP solution (dense, size numCols), used to score cuts.
 * @param numCols number of columns of the model `x` refers to.
 * @param generatorTag short name of the calling generator (e.g. "Gomory"),
 *        used only for CBC_CLIQUE_POOL_DEBUG-style diagnostic output.
 * @return true if any cuts were removed.
 **/
bool cbcFilterGeneratedCuts(OsiCuts &cs, int firstRowCut, const double *x,
  int numCols, const char *generatorTag);

#endif //CBCCUTPOOLFILTER_HPP
