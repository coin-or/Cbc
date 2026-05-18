/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for extending cliques
 *
 * @file CoinCliqueExtender.cpp
 * @brief Clique extender
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>
#include "CoinCliqueExtender.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinCliqueList.hpp"
#include "CoinStaticConflictGraph.hpp"

#define CLQEXT_EPS 1e-6

static void shuffle_array(size_t *arr, size_t n);

struct CompareCost {
    explicit CompareCost(const double *costs) { this->costs_ = costs; }

    bool operator () (const size_t &i, const size_t &j) {
        if(fabs(costs_[i] - costs_[j]) >= CLQEXT_EPS) {
            return costs_[i] + CLQEXT_EPS <= costs_[j];
        }

        return i < j;
    }

    const double *costs_;
};

CoinCliqueExtender::CoinCliqueExtender(const CoinConflictGraph *cgraph, size_t extMethod, const double *rc, const double maxRC)
: maxCandidates_(512)
{
    const size_t cgSize = cgraph->size();
    cgraph_ = cgraph;
    extendedCliques_ = new CoinCliqueList(4096, 32768);
    nCandidates_ = nNewClique_ = 0;
    candidates_ = std::vector<size_t>(cgSize);
    newClique_ = std::vector<size_t>(cgSize);

    iv_ = std::vector<char>(cgSize);
    iv2_ = std::vector<char>(cgSize);

    extMethod_ = extMethod;
    rc_ = rc;
    maxRC_ = maxRC;

    if ((extMethod_ == 4 || extMethod_ == 5) && !rc_) {
        fprintf(stderr, "Warning: using random selection for extension since no costs were informed.\n");
        extMethod_ = 1;
    }

    switch (extMethod) {
        case 0: //no extension
            //nothing to do
            break;
        case 1: { //random extension
            //nothing to do
            break;
        }
        case 2: { //max degree extension
            costs_ = std::vector<double>(cgSize);
            for (size_t i = 0; i < cgSize; i++) {
#ifdef DEBUGCG
                assert(cgraph_->degree(i) < cgSize);
#endif
                costs_[i] = ((double) (cgSize - cgraph_->degree(i)));
            }
        }
            break;
        case 3: {//modified degree extension
            double maxModDegree = 0.0;
            costs_ = std::vector<double>(cgSize);
            for (size_t i = 0; i < cgSize; i++) {
                costs_[i] = cgraph_->modifiedDegree(i);
                maxModDegree = std::max(maxModDegree, costs_[i]);
            }
            for (size_t i = 0; i < cgSize; i++) {
                costs_[i] = maxModDegree - costs_[i] + 1.0;
#ifdef DEBUGCG
                assert(costs_[i] >= 1.0 && costs_[i] <= maxModDegree + 1.0);
#endif
            }
            break;
        }
        case 4: { //priority greedy extension
            //nothing to do - costs are the reduced cost array
            break;
        }
        case 5: { //reduced cost + modified degree
            costs_ = std::vector<double>(cgSize);
            double minRCost, maxRCost;
            size_t minMDegree, maxMDegree;
            minRCost = maxRCost = rc[0];
            minMDegree = maxMDegree = cgraph_->modifiedDegree(0);
            for (size_t i = 1; i < cgSize; i++) {
                minRCost = std::min(minRCost, rc_[i]);
                maxRCost = std::max(maxRCost, rc_[i]);
                minMDegree = std::min(minMDegree, cgraph_->modifiedDegree(i));
                maxMDegree = std::max(maxMDegree, cgraph_->modifiedDegree(i));
            }
            for (size_t i = 0; i < cgSize; i++) {
                const double normRC = (rc_[i] - minRCost) / (maxRCost - minRCost + CLQEXT_EPS);
                const double normMD = 1.0 - (((double)(cgraph_->modifiedDegree(i) - minMDegree))
                                    / ((double)(maxMDegree - minMDegree) + CLQEXT_EPS));
                costs_[i] = (700.0 * normRC) + (300.0 * normMD);
#ifdef DEBUGCG
                assert(normRC >= 0.0 && normRC <= 1.0);
                assert(normMD >= 0.0 && normMD <= 1.0);
                assert(costs_[i] >= 0.0 && costs_[i] <= 1000.0);
#endif
            }
            break;
        }
        default:
            fprintf(stderr, "Invalid option %lu\n", extMethod_);
            abort();
    }
}

CoinCliqueExtender::~CoinCliqueExtender() {
    delete extendedCliques_;
}

bool CoinCliqueExtender::randomExtension(const size_t *clqIdxs, const size_t clqSize) {
#ifdef DEBUGCG
    assert(clqSize > 0);
#endif

    fillCandidates(clqIdxs, clqSize);

    if (nCandidates_ == 0) {
        return false;
    }

    shuffle_array(candidates_.data(), nCandidates_);
    nCandidates_ = std::min(nCandidates_, maxCandidates_);

    for (size_t i = 0; i < nCandidates_; i++) {
        /* need to have conflict with all nodes in clique */
        const size_t selected = candidates_[i];
        bool insert = true;
        for (size_t j = clqSize; j < nNewClique_; j++) {
            if (!cgraph_->conflicting(newClique_[j], selected)) {
                insert = false;
                break;
            }
        }
        if (insert) {
            newClique_[nNewClique_++] = selected;
        }
    }

    if (nNewClique_ == clqSize) {
        return false;
    }

#ifdef DEBUGCG
    CoinCliqueList::validateClique(cgraph_, newClique_.data(), nNewClique_);
#endif

    extendedCliques_->addClique(nNewClique_, newClique_.data());

    return true;
}

bool CoinCliqueExtender::greedySelection(const size_t *clqIdxs, const size_t clqSize, const double *costs) {
#ifdef DEBUGCG
    assert(clqSize > 0);
#endif

    fillCandidates(clqIdxs, clqSize);

    if (nCandidates_ == 0) {
        return false;
    }

    const size_t n = nCandidates_;
    nCandidates_ = std::min(nCandidates_, maxCandidates_);
    std::partial_sort(candidates_.begin(), candidates_.begin() + nCandidates_, candidates_.begin() + n, CompareCost(costs));

    for (size_t i = 0; i < nCandidates_; i++) {
        /* need to have conflict with all nodes in clique */
        const size_t selected = candidates_[i];
        bool insert = true;
        for (size_t j = clqSize; j < nNewClique_; j++) {
            if (!cgraph_->conflicting(newClique_[j], selected)) {
                insert = false;
                break;
            }
        }
        if (insert) {
            newClique_[nNewClique_++] = selected;
        }
    }

    if (nNewClique_ == clqSize) {
        return false;
    }

#ifdef DEBUGCG
    CoinCliqueList::validateClique(cgraph_, newClique_.data(), nNewClique_);
#endif

    extendedCliques_->addClique(nNewClique_, newClique_.data());

    return true;
}

bool compareCliqueWeight(const std::pair<size_t, double> &e1, const std::pair<size_t, double> &e2) {
    if (fabs(e1.second - e2.second) >= CLQEXT_EPS) {
        return e1.second >= e2.second + CLQEXT_EPS;
    }

    return e1.first < e2.first;
}

bool CoinCliqueExtender::extendClique(const size_t *clqIdxs, const size_t clqSize) {
    if (extMethod_ == 0) {
        return false;
    }

#ifdef DEBUGCG
    if(extMethod_ == 4 || extMethod_ == 5) {
        assert(rc_);
    }
#endif

    bool result = false;

    switch (extMethod_) {
        case 1: //random
            result = randomExtension(clqIdxs, clqSize);
            break;
        case 2: //max degree
            result = greedySelection(clqIdxs, clqSize, costs_.data());
            break;
        case 3: //modified degree
            result = greedySelection(clqIdxs, clqSize, costs_.data());
            break;
        case 4: //priority greedy (reduced cost)
            result = greedySelection(clqIdxs, clqSize, rc_);
            break;
        case 5: //reduced cost + modified degree
            result = greedySelection(clqIdxs, clqSize, costs_.data());
            break;
        default:
            fprintf(stderr, "Invalid option %lu\n", extMethod_);
            abort();
    }

    return result;
}

void CoinCliqueExtender::fillCandidates(const size_t *clqIdxs, const size_t clqSize) {
    const size_t cgSize = cgraph_->size();
    size_t nodeSD = 0, minDegree = cgSize;

    nNewClique_ = nCandidates_ = 0;

    /* picking node with the smallest degree */
    /* adding clique elements in newClique */
    for (size_t i = 0; i < clqSize; i++) {
        const size_t clqIdx = clqIdxs[i];
        const size_t degree = cgraph_->degree(clqIdx);

        if (degree < minDegree) {
            minDegree = degree;
            nodeSD = clqIdx;
        }

        newClique_[nNewClique_++] = clqIdx;
        iv_[clqIdx] = 1;
    }

    const std::pair<size_t, const size_t*> rescg = cgraph_->conflictingNodes(nodeSD, candidates_.data(), iv2_.data());

    for (size_t i = 0; i < rescg.first; i++) {
        const size_t node = rescg.second[i];

        if (iv_[node]) { //already inserted at clique
            continue;
        }

        if (rc_ && rc_[node] >= maxRC_ + CLQEXT_EPS) {
            continue;
        }

        bool insert = true;
        for (size_t j = 0; j < clqSize; j++) {
            if (!cgraph_->conflicting(node, clqIdxs[j])) {
                insert = false;
                break;
            }
        }
        if (insert) {
            candidates_[nCandidates_++] = node;
        }
    }

    for (size_t i = 0; i < clqSize; i++) {
        iv_[clqIdxs[i]] = 0;
    }

#ifdef DEBUGCG
    assert(nCandidates_ <= rescg.first);
#endif
}

void CoinCliqueExtender::setMaxCandidates(const size_t maxCandidates) {
    maxCandidates_ = maxCandidates;
}

size_t CoinCliqueExtender::nCliques() const {
    if (!extendedCliques_) {
        return 0;
    }
    return extendedCliques_->nCliques();
}

const size_t* CoinCliqueExtender::getClique(const size_t i) const {
#ifdef DEBUGCG
    assert(i < extendedCliques_->nCliques());
#endif
    return extendedCliques_->cliqueElements(i);
}

size_t CoinCliqueExtender::getCliqueSize(const size_t i) const {
#ifdef DEBUGCG
    assert(i < extendedCliques_->nCliques());
#endif
    return extendedCliques_->cliqueSize(i);
}

static void shuffle_array (size_t *arr, size_t n) {
    if (n <= 1) {
        return;
    }

    size_t tmp;

    srand (time(NULL));

    for (size_t i = n - 1; i > 0; i--) {
        size_t j = rand() % (i + 1);
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}
