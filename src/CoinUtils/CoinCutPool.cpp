/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for storing a pool of cuts, removing
 * the repeated and dominated ones. It also filters the cuts
 * according to their scores. Given a solution x* of the LP
 * relaxation of a MILP, the score of a cut C is calculated as
 * S(C) = viol(C) / actv(C), where viol(C) is the violation
 * of the cut with respect to x* and actv(C) is the
 * number of variables in C whose values x* are greater than
 * zero. An auxiliary array is used to identify, for each variable,
 * the cut in the pool with the best score that contains this variable.
 * A cut is only inserted into the cut pool if it has the best score for
 * at least one variable.
 *
 * @file CoinCutPool.cpp
 * @brief Class for storing a pool of cuts
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include "CoinCutPool.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <algorithm>

#define CUTPOOL_EPS 1e-8

struct CompareIdxs {
    explicit CompareIdxs(const int *idxs) { this->idxs_ = idxs; }

    bool operator () (const size_t &i, const size_t &j) {
        return idxs_[i] < idxs_[j];
    }

    const int *idxs_;
};

CoinCut::CoinCut(const int *idxs, const double *coefs, size_t nz, double rhs)
  : idxs_(std::vector<int>(nz))
  , coefs_(std::vector<double>(nz))
  , rhs_(rhs)
{
    for (size_t i = 0; i < nz; i++) {
        idxs_[i] = i;
    }

    std::sort(idxs_.begin(), idxs_.end(), CompareIdxs(idxs));

    for (size_t i = 0; i < nz; i++) {
        const int pos = idxs_[i];
        idxs_[i] = idxs[pos];
        coefs_[i] = coefs[pos];
    }

#ifdef DEBUGCG
    assert(std::is_sorted(idxs_.begin(), idxs_.end()));
#endif
}

CoinCut::~CoinCut() {}

int binarySearch(const int *v, const int n, const int x) {
    int mid, left = 0, right = n - 1;

    while (left <= right) {
        mid = (left + right) / 2;

        if (v[mid] == x) {
            return mid;
        } else if (v[mid] < x) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return -1;
}

bool CoinCut::dominates(const CoinCut *other) const {
    /* thisRHS == 0 && otherRHS < 0 */
    if (fabs(this->rhs()) <= CUTPOOL_EPS && other->rhs() <= -CUTPOOL_EPS) {
        return false;
    }

    /* thisRHS > 0 && otherRHS == 0 */
    if (this->rhs() >= CUTPOOL_EPS && fabs(other->rhs()) <= CUTPOOL_EPS) {
        return false;
    }

    int sizeA = this->size(), sizeB = other->size();
    const int *idxsA = this->idxs(), *idxsB = other->idxs();
    const double *coefsA = this->coefs(), *coefsB = other->coefs();
    const double rhsA = this->rhs(), rhsB = other->rhs();
    double normConstA, normConstB;

    std::vector<bool> iv(sizeB, false);

    if (fabs(rhsA) <= CUTPOOL_EPS) {
        normConstA = 1.0;
    } else {
        normConstA = rhsA;
    }

    if (fabs(rhsB) <= CUTPOOL_EPS) {
        normConstB = 1.0;
    } else {
        normConstB = rhsB;
    }

    for (int i = 0; i < sizeA; i++) {
        const int idxA = idxsA[i];
        const double normCoefA = (coefsA[i] / normConstA) * 1000.0;
        int posB = binarySearch(idxsB, sizeB, idxA);

        if (posB == -1 && normCoefA <= -CUTPOOL_EPS) {
            return false;
        } else if (posB != -1) {
            double normCoefB = (coefsB[posB] / normConstB) * 1000.0;
            iv[posB] = true;
            if (normCoefA + CUTPOOL_EPS <= normCoefB) {
                return false;
            }
        }
    }

    for (int i = 0; i < sizeB; i++) {
        if (iv[i]) {
            continue;
        }

        int idxB = idxsB[i];
        double normCoefB = (coefsB[i] / normConstB) * 1000.0;
        int posA = binarySearch(idxsA, sizeA, idxB);

        if (posA == -1 && normCoefB >= CUTPOOL_EPS) {
            return false;
        }
    }

    return true;
}

CoinCutPool::CoinCutPool(const double *x, int numCols) {
    x_ = x;
    nullCuts_ = 0;

    bestCutByCol_ = std::vector<int>(numCols, -1);

    nCuts_ = 0;
    cutsCap_ = 1024;
    cuts_ = std::vector<CoinCut*>(cutsCap_);
    cutFrequency_ = std::vector<size_t>(cutsCap_);
    cutFitness_ = std::vector<double>(cutsCap_);
}

CoinCutPool::~CoinCutPool() {
    for (size_t i = 0; i < nCuts_; i++) {
        delete cuts_[i];
    }
}

bool CoinCutPool::add(const int *idxs, const double *coefs, int nz, double rhs) {
    CoinCut *cut = new CoinCut(idxs, coefs, nz, rhs);

    checkMemory();

    if (updateCutFrequency(cut) == 0) {
        delete cut;
        return false;
    }

    cuts_[nCuts_++] = cut;
    return true;
}

size_t CoinCutPool::updateCutFrequency(const CoinCut *cut) {
    const int nz = cut->size();
    const int *idxs = cut->idxs();
    const double fitness = calculateFitness(cut);

    cutFitness_[nCuts_] = fitness;
    cutFrequency_[nCuts_] = 0;

    for (int i = 0; i < nz; i++) {
        const int idx = idxs[i];

        if (bestCutByCol_[idx] == -1) {
            bestCutByCol_[idx] = nCuts_;
            cutFrequency_[nCuts_]++;
        } else {
            const int currBestCut = bestCutByCol_[idx];
            const double currFitness = cutFitness_[currBestCut];

            if (fitness >= currFitness + CUTPOOL_EPS) {
                cutFrequency_[nCuts_]++;
                bestCutByCol_[idx] = nCuts_;
#ifdef DEBUGCG
                assert(cutFrequency_[currBestCut] > 0);
#endif
                cutFrequency_[currBestCut]--;

                if (cutFrequency_[currBestCut] == 0) {
                    delete cuts_[currBestCut];
                    cuts_[currBestCut] = NULL;
                    nullCuts_++;
                }
            }
        }
    }

    return cutFrequency_[nCuts_];
}

double CoinCutPool::calculateFitness(const CoinCut *cut) const {
    const int size = cut->size();
    const int *idxs = cut->idxs();
    const double *coefs = cut->coefs();
    const double rhs = cut->rhs();
    double lhs = 0.0, minCoef, maxCoef;
    size_t activeCols = 0;

    minCoef = maxCoef = coefs[0];

    for (size_t i = 0; i < size; i++) {
#ifdef DEBUGCG
        assert(fabs(coefs[i]) >= CUTPOOL_EPS);
#endif

        if (fabs(x_[idxs[i]]) >= CUTPOOL_EPS) {
            lhs += (coefs[i] * x_[idxs[i]]);
            activeCols++;
            minCoef = std::min(minCoef, coefs[i]);
            maxCoef = std::max(maxCoef, coefs[i]);
        }
    }

    const double diffOfCoefs = fabs(maxCoef - minCoef) + fabs(maxCoef - rhs) + fabs(minCoef - rhs);
    const double violation = lhs - rhs;

#ifdef DEBUGCG
    assert(violation >= CUTPOOL_EPS);
    assert(activeCols > 0);
#endif

    return ((violation / ((double) activeCols)) * 100000.0) + ((1.0 / (diffOfCoefs + 1.0)) * 1000.0);
}

void CoinCutPool::checkMemory() {
    if (nCuts_ + 1 <= cutsCap_) {
        return;
    }

    cutsCap_ *= 2;
    cuts_.resize(cutsCap_);
    cutFrequency_.resize(cutsCap_);
    cutFitness_.resize(cutsCap_);
}

void CoinCutPool::removeDominated() {
    if (nCuts_ <= 1) {
        return;
    }

    if (nullCuts_) {
        removeNullCuts();
    }

#ifdef DEBUGCG
    assert(nullCuts_ == 0);
#endif

    for (size_t i = 0; i < nCuts_; i++) {
#ifdef DEBUGCG
        assert(cuts_[i]);
#endif
        for (size_t j = i + 1; j < nCuts_; j++) {
#ifdef DEBUGCG
            assert(cuts_[j]);
#endif
            int chkDm = checkCutDomination(i, j);

            if (chkDm == 0) { //cut i dominates cut j
                delete cuts_[j];
                cuts_[j] = NULL;
                nullCuts_++;
                continue;
            } else if (chkDm == 1) { //cut j dominates cut i
                delete cuts_[i];
                cuts_[i] = NULL;
                nullCuts_++;
                break;
            }
        }
    }

#ifdef DEBUGCG
    assert(nullCuts_ >= 0 && nullCuts_ < nCuts_);
#endif

    if(nullCuts_ > 0) {
        removeNullCuts();
    }
}

int CoinCutPool::checkCutDomination(size_t idxA, size_t idxB) {
    const CoinCut *cutA = cuts_[idxA];
    const CoinCut *cutB = cuts_[idxB];

    /* checks if cutA dominates cutB */
    if (cutA->dominates(cutB)) {
        return 0;
    }

    /* checks if cutB dominates cutA */
    if (cutB->dominates(cutA)) {
        return 1;
    }

    /* cutA and cutB are not dominated */
    return 2;
}

size_t CoinCutPool::numCuts() const {
    return nCuts_;
}

const int* CoinCutPool::cutIdxs(size_t i) const {
#ifdef DEBUGCG
    assert(i < nCuts_);
    assert(cuts_[i]);
#endif
    return cuts_[i]->idxs();
}

const double* CoinCutPool::cutCoefs(size_t i) const {
#ifdef DEBUGCG
    assert(i < nCuts_);
    assert(cuts_[i]);
#endif
    return cuts_[i]->coefs();
}

int CoinCutPool::cutSize(size_t i) const {
#ifdef DEBUGCG
    assert(i < nCuts_);
    assert(cuts_[i]);
#endif
    return cuts_[i]->size();
}

double CoinCutPool::cutRHS(size_t i) const {
#ifdef DEBUGCG
    assert(i < nCuts_);
    assert(cuts_[i]);
#endif
    return cuts_[i]->rhs();
}

void CoinCutPool::removeNullCuts() {
    if (nCuts_ == 0) {
        return;
    }

    if (nullCuts_ > 0) {
        size_t last = 0;

        for (size_t i = 0; i < nCuts_; i++) {
            if (cuts_[i]) {
                cuts_[last++] = cuts_[i];
            }
        }

#ifdef DEBUGCG
        assert(last == (nCuts_ - nullCuts_));
#endif

        nCuts_ -= nullCuts_;
        nullCuts_ = 0;
    }
}

