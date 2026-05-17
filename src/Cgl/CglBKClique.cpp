/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for separating clique cuts.
 * It uses the Bron-Kerbosch (BK) algorithm.
 *
 * @file CglBKClique.hpp
 * @brief Clique cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <limits>
#include <utility>
#include <vector>
#include <OsiCuts.hpp>
#include <OsiRowCut.hpp>
#include <CoinTime.hpp>

#include "CglBKClique.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinCliqueList.hpp"
#include "CoinCliqueExtender.hpp"
#include "CoinCutPool.hpp"

#define BKCLQ_MULTIPLIER 1000.0
#define BKCLQ_EPS 1e-6

std::atomic< size_t > CglBKClique::sepCuts_(0);
std::atomic< double > CglBKClique::sepTime_(0.0);

static void *xmalloc(const size_t size);

CglBKClique::CglBKClique()
  : cap_(0)
  , minFrac_(0.001)
  , minViol_(0.02)
  , pivotingStrategy_(CoinBronKerbosch::PivotingStrategy::Weight)
  , extMethod_(4)
  , maxCallsBK_(1000)
{
  minWeight_ = floor(BKCLQ_MULTIPLIER + (minViol_ * BKCLQ_MULTIPLIER));
  vertexWeight_ = NULL;
  rc_ = NULL;
  idxs_ = NULL;
  idxMap_ = NULL;
  coefs_ = NULL;
  inducedVert_ = NULL;
  currClq_ = NULL;
  callsBK_ = 0;
  completeBK_ = false;
  maxSeconds_ = 0.0;
  maxInducedSize_ = 10000;
}

CglBKClique::CglBKClique(const CglBKClique &rhs)
  : CglCutGenerator(rhs)
{
  this->cap_ = rhs.cap_;

  if (this->cap_ > 0) {
    this->vertexWeight_ = (double *)xmalloc(sizeof(double) * this->cap_ * 2);
    this->rc_ = (double *)xmalloc(sizeof(double) * this->cap_ * 2);
    this->idxs_ = (int *)xmalloc(sizeof(int) * this->cap_);
    this->idxMap_ = (int *)xmalloc(sizeof(int) * this->cap_);
    this->coefs_ = (double *)xmalloc(sizeof(double) * this->cap_);
    this->inducedVert_ = (size_t *)xmalloc(sizeof(size_t) * this->cap_ * 2);
    this->currClq_ = (size_t *)xmalloc(sizeof(size_t) * this->cap_ * 2);
  } else {
    this->vertexWeight_ = NULL;
    this->rc_ = NULL;
    this->idxs_ = NULL;
    this->idxMap_ = NULL;
    this->coefs_ = NULL;
    this->inducedVert_ = NULL;
    this->currClq_ = NULL;
  }

  this->maxCallsBK_ = rhs.maxCallsBK_;
  this->extMethod_ = rhs.extMethod_;
  this->minFrac_ = rhs.minFrac_;
  this->minViol_ = rhs.minViol_;
  this->pivotingStrategy_ = rhs.pivotingStrategy_;
  this->minWeight_ = rhs.minWeight_;
  this->callsBK_ = rhs.callsBK_;
  this->completeBK_ = rhs.completeBK_;
  this->maxSeconds_ = rhs.maxSeconds_;
  this->maxInducedSize_ = rhs.maxInducedSize_;
}

CglBKClique::~CglBKClique()
{
  if (this->vertexWeight_) {
    free(this->vertexWeight_);
  }
  if (this->rc_) {
    free(this->rc_);
  }
  if (this->idxs_) {
    free(this->idxs_);
  }
  if (this->idxMap_) {
    free(this->idxMap_);
  }
  if (this->coefs_) {
    free(this->coefs_);
  }
  if (this->inducedVert_) {
    free(this->inducedVert_);
  }
  if (this->currClq_) {
    free(this->currClq_);
  }
}

void CglBKClique::refreshSolver(OsiSolverInterface *solver)
{
  solver->checkCGraph();
  // Get integer information
  solver->getColType(true);
}

CglCutGenerator *CglBKClique::clone() const
{
  return new CglBKClique(*this);
}

void CglBKClique::generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info)
{
  if (si.getNumCols() == 0 || si.getNumRows() == 0) {
    return;
  }

  double startSep = CoinCpuTime();
  const CoinConflictGraph *cgraph = si.getCGraph();

  if (si.getNumCols() != cgraph->size() / 2) {
    fprintf(stderr, "Invalid conflict graph! Number of columns %d ... in graph %lu\n",
      si.getNumCols(), cgraph->size() / 2);
    abort();
  }

  checkMemory(si.getNumCols());

  // When a wall-clock deadline is set, adaptively reduce the BK call limit
  // if time is already tight before we start.  This prevents the recursive
  // BK algorithm from running for hundreds of seconds on dense graphs.
  const double wallDeadline = (maxSeconds_ > 0.0)
    ? CoinGetTimeOfDay() + maxSeconds_
    : std::numeric_limits< double >::max();
  if (maxSeconds_ > 0.0 && maxSeconds_ < 30.0) {
    // Very little time left: drastically limit BK recursion depth.
    const size_t cappedCalls = std::max(size_t(10),
      size_t(maxCallsBK_ * maxSeconds_ / 30.0));
    maxCallsBK_ = std::min(maxCallsBK_, cappedCalls);
  }

  CoinCliqueList *initialCliques = separateCliques(si);

  if (initialCliques->nCliques() > 0) {
    if (!extMethod_ || CoinGetTimeOfDay() >= wallDeadline) {
      insertCuts(si, info, initialCliques, cs);
    } else {
      CoinCliqueList *extCliques = extendCliques(si, initialCliques);
      insertCuts(si, info, extCliques, cs);
      delete extCliques;
    }
  }

  delete initialCliques;
  {
    const double delta = CoinCpuTime() - startSep;
    double old = CglBKClique::sepTime_.load(std::memory_order_relaxed);
    while (!CglBKClique::sepTime_.compare_exchange_weak(old, old + delta, std::memory_order_relaxed))
      ;
  }
}

void CglBKClique::checkMemory(const size_t newNumCols)
{
  if (cap_ < newNumCols) {
    if (cap_ > 0) {
#ifdef DEBUGCG
      assert(vertexWeight_);
      assert(rc_);
      assert(idxs_);
      assert(idxMap_);
      assert(coefs_);
      assert(inducedVert_);
      assert(currClq_);
#endif
      free(vertexWeight_);
      free(rc_);
      free(idxs_);
      free(idxMap_);
      free(coefs_);
      free(inducedVert_);
      free(currClq_);
    }

    vertexWeight_ = (double *)xmalloc(sizeof(double) * newNumCols * 2);
    rc_ = (double *)xmalloc(sizeof(double) * newNumCols * 2);
    idxs_ = (int *)xmalloc(sizeof(int) * newNumCols);
    idxMap_ = (int *)xmalloc(sizeof(int) * newNumCols);
    coefs_ = (double *)xmalloc(sizeof(double) * newNumCols);
    inducedVert_ = (size_t *)xmalloc(sizeof(size_t) * newNumCols * 2);
    currClq_ = (size_t *)xmalloc(sizeof(size_t) * newNumCols * 2);
    cap_ = newNumCols;
  }
}

CoinCliqueList *CglBKClique::separateCliques(const OsiSolverInterface &si)
{
  const size_t numCols = si.getNumCols();
  const double *x = si.getColSolution();
  CoinCliqueList *initialCliques = new CoinCliqueList(4096, 32768);
  const CoinConflictGraph *cgraph = si.getCGraph();
  size_t n = 0;

  // generating the subgraph induced by the fractional variables and variables at one
  for (size_t i = 0; i < numCols; i++) { // variables
    const size_t degree = cgraph->degree(i);
    if (degree < 2) {
      // disconsidering variables that have no conflicts or just conflicts involving their complements
      continue;
    } else if (x[i] + BKCLQ_EPS <= minFrac_) {
      // variables at zero in x are disconsidered
      continue;
    } else {
      inducedVert_[n] = i;
      vertexWeight_[n] = x[i] * BKCLQ_MULTIPLIER;
      n++;
    }
  }
  for (size_t i = numCols; i < numCols * 2; i++) { // complements
    const size_t degree = cgraph->degree(i);
    const double xc = 1.0 - x[i - numCols];
    if (degree < 2) {
      // disconsidering variables that have no conflicts or just conflicts involving their complements
      continue;
    } else if (xc + BKCLQ_EPS <= minFrac_) {
      // variables at zero in x are disconsidered
      continue;
    } else {
      inducedVert_[n] = i;
      vertexWeight_[n] = xc * BKCLQ_MULTIPLIER;
      n++;
    }
  }

  if (n == 0) {
    return initialCliques;
  }

  // Cap the induced subgraph to avoid O(n²) cost in CoinBronKerbosch when
  // many fractional variables (originals + complements) are present.
  // Vertices are ranked by LP-value × sqrt(conflict-degree): this jointly
  // favours high-LP variables (efficient violators of small cliques) and
  // high-degree variables (members of large cliques).
  // Variables near x = 0.5 appear twice (original weight ≈ 0.5 and
  // complement weight ≈ 0.5) and both receive the same moderate priority,
  // so neither form is unfairly dropped.
  if (maxInducedSize_ > 0 && n > maxInducedSize_) {
    std::vector< std::pair< double, size_t > > scored(n);
    for (size_t i = 0; i < n; i++) {
      const double deg = static_cast< double >(cgraph->degree(inducedVert_[i]) + 1);
      scored[i] = { vertexWeight_[i] * std::sqrt(deg), i };
    }
    std::partial_sort(scored.begin(), scored.begin() + maxInducedSize_,
      scored.end(),
      [](const std::pair< double, size_t > &a,
        const std::pair< double, size_t > &b) { return a.first > b.first; });
    std::vector< size_t > tmpVert(maxInducedSize_);
    std::vector< double > tmpWeight(maxInducedSize_);
    for (size_t i = 0; i < maxInducedSize_; i++) {
      const size_t src = scored[i].second;
      tmpVert[i] = inducedVert_[src];
      tmpWeight[i] = vertexWeight_[src];
    }
    std::copy(tmpVert.begin(), tmpVert.end(), inducedVert_);
    std::copy(tmpWeight.begin(), tmpWeight.end(), vertexWeight_);
    n = maxInducedSize_;
  }

  CoinConflictGraph *ppcg = new CoinStaticConflictGraph(cgraph, n, inducedVert_);
#ifdef DEBUGCG
  assert(ppcg->size() == n);
#endif

  completeBK_ = true;

  if (ppcg->size() >= 2) {
    ppcg->computeModifiedDegree();
    CoinBronKerbosch *bk = new CoinBronKerbosch(ppcg, vertexWeight_, pivotingStrategy_);
    bk->setMaxCalls(maxCallsBK_);
    bk->setMinWeight(minWeight_);
    bk->findCliques();
    completeBK_ = bk->completedSearch();
    maxCallsBK_ = bk->numCalls();

    if (bk->nCliques() > 0) {
      for (size_t i = 0; i < bk->nCliques(); i++) {
        const size_t cliqueSize = bk->getCliqueSize(i);
        const size_t *cliqueIdxs = bk->getClique(i);
        for (size_t j = 0; j < cliqueSize; j++) {
          currClq_[j] = inducedVert_[cliqueIdxs[j]];
        }
        initialCliques->addClique(cliqueSize, currClq_);
      }
    }

    delete bk;
  }

  delete ppcg;

  return initialCliques;
}

CoinCliqueList *CglBKClique::extendCliques(const OsiSolverInterface &si, const CoinCliqueList *initialCliques)
{
  CoinCliqueList *extCliques = new CoinCliqueList(4096, 32768);
  const double *rCost = si.getReducedCost();
  const size_t numCols = si.getNumCols();
  const CoinConflictGraph *cgraph = si.getCGraph();

  // setting reduced costs
  for (size_t i = 0; i < numCols; i++) {
    rc_[i] = rCost[i];
    rc_[i + numCols] = -rc_[i];
  }

  CoinCliqueExtender clqe(cgraph, extMethod_, rc_, 100.0);

  for (size_t i = 0; i < initialCliques->nCliques(); i++) {
    const size_t *clqEl = initialCliques->cliqueElements(i);
    const size_t nClqEl = initialCliques->cliqueSize(i);
    const bool extended = clqe.extendClique(clqEl, nClqEl);

    /* adds clique if it is not extended */
    if (!extended) {
      extCliques->addClique(nClqEl, clqEl);
    }
  }

  /* adding all extended cliques */
  for (size_t i = 0; i < clqe.nCliques(); i++) {
    extCliques->addClique(clqe.getCliqueSize(i), clqe.getClique(i));
  }

  return extCliques;
}

void CglBKClique::insertCuts(const OsiSolverInterface &si, const CglTreeInfo &info, const CoinCliqueList *cliques, OsiCuts &cs)
{
  const double *x = si.getColSolution();
  const size_t numCols = si.getNumCols();
  CoinCutPool cutpool(x, numCols);

  for (size_t i = 0; i < cliques->nCliques(); i++) {
    const size_t clqSize = cliques->cliqueSize(i);
    const size_t *el = cliques->cliqueElements(i);
    double rhs = 1.0;
    int cutSize = 0;
    size_t duplicated = 0;

    std::fill(idxMap_, idxMap_ + numCols, -1);

    for (size_t j = 0; j < clqSize; j++) {
      if (el[j] < numCols) {
        if (idxMap_[el[j]] == -1) {
          idxMap_[el[j]] = cutSize;
          idxs_[cutSize] = ((int)el[j]);
          coefs_[cutSize] = 1.0;
          cutSize++;
        } else {
          coefs_[idxMap_[el[j]]] += 1.0;
          assert(coefs_[idxMap_[el[j]]] == 0.0);
          duplicated++;
        }
      } else {
        rhs -= 1.0;
        if (idxMap_[el[j] - numCols] == -1) {
          idxMap_[el[j] - numCols] = cutSize;
          idxs_[cutSize] = ((int)(el[j] - numCols));
          coefs_[cutSize] = -1.0;
          cutSize++;
        } else {
          coefs_[idxMap_[el[j] - numCols]] -= 1.0;
          assert(coefs_[idxMap_[el[j] - numCols]] == 0.0);
          duplicated++;
        }
      }
    }

    assert(duplicated == 0 || duplicated == 1);
    if (duplicated == 1) {
      int last = 0;
      rhs = 0.0;
      for (int k = 0; k < cutSize; k++) {
        assert(coefs_[k] == -1.0 || coefs_[k] == 0.0 || coefs_[k] == 1.0);
        if (coefs_[k] == -1.0 || coefs_[k] == 1.0) {
          idxs_[last] = idxs_[k];
          coefs_[last] = coefs_[k];
          last++;
          if (coefs_[k] == -1.0) {
            rhs -= 1.0;
          }
        }
      }
      cutSize = last;
    }

    cutpool.add(idxs_, coefs_, cutSize, rhs);
  }

  cutpool.removeNullCuts();

  const size_t numberRowCutsBefore = cs.sizeRowCuts();
  for (size_t i = 0; i < cutpool.numCuts(); i++) {
    osrc_.setRow(cutpool.cutSize(i), cutpool.cutIdxs(i), cutpool.cutCoefs(i));
    osrc_.setUb(cutpool.cutRHS(i));
    cs.insertIfNotDuplicate(osrc_);
  }

  size_t numberRowCutsAfter = cs.sizeRowCuts();
  CglBKClique::sepCuts_.fetch_add(numberRowCutsAfter - numberRowCutsBefore, std::memory_order_relaxed);

  if (!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
    numberRowCutsAfter = cs.sizeRowCuts();
    for (size_t i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
      cs.rowCutPtr(i)->setGloballyValid();
    }
  }
}

void CglBKClique::setMaxCallsBK(size_t maxCallsBK)
{
  this->maxCallsBK_ = maxCallsBK;
}

void CglBKClique::setExtendingMethod(size_t extMethod)
{
  if (extMethod > 6) {
    fprintf(stderr, "Invalid value for parameter extMethod (%ld).\n", (size_t)extMethod);
    abort();
  }

  this->extMethod_ = extMethod;
}

void CglBKClique::setMinFrac(const double minFrac)
{
  minFrac_ = minFrac;
}

void CglBKClique::setMinViol(const double minViol)
{
  minViol_ = minViol;
  minWeight_ = BKCLQ_MULTIPLIER + (minViol_ * BKCLQ_MULTIPLIER);
}

void CglBKClique::setPivotingStrategy(const CoinBronKerbosch::PivotingStrategy pivotingStrategy)
{
  pivotingStrategy_ = pivotingStrategy;
}

static void *xmalloc(const size_t size)
{
  void *result = malloc(size);
  if (!result) {
    fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
    abort();
  }

  return result;
}
