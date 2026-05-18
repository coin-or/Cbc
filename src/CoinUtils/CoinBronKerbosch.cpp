/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for finding cliques with weights
 * greater than a certain threshold. It implements the
 * Bron-Kerbosch Algorithm.
 *
 * @file CoinBronKerbosch.cpp
 * @brief Bron-Kerbosch Algorithm
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "CoinBronKerbosch.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinCliqueList.hpp"

#define INT_SIZE ((size_t)(8 * sizeof(int)))
#define BK_EPS 1e-6

static void shuffle_vertices (BKVertex *vertices, size_t n);

bool compareNodes(const BKVertex &u, const BKVertex &v) {
    return u.fitness >= v.fitness + BK_EPS;
}

CoinBronKerbosch::CoinBronKerbosch(const CoinConflictGraph *cgraph, const double *weights, PivotingStrategy pivotingStrategy)
  : minWeight_(0.0)
  , calls_(0)
  , maxCalls_(std::numeric_limits<size_t >::max())
  , pivotingStrategy_(pivotingStrategy)
  , completeSearch_(false)
{
    const size_t cgSize = cgraph->size();
    size_t maxDegree = 0;

    cliques_ = NULL;
    cgraph_ = cgraph;
    vertices_ = std::vector<BKVertex>();
    vertices_.reserve(cgSize);

    //filling information about vertices
    for (size_t u = 0; u < cgSize; u++) {
        const size_t degree = cgraph_->degree(u);

        if (degree > 0) {
            BKVertex bkv;
            bkv.idx = u;
            bkv.weight = weights[u];
            bkv.degree = degree;
            vertices_.push_back(bkv);
            maxDegree = std::max(maxDegree, degree);
        }
    }

    if (vertices_.empty()) {
        return;
    }

    computeFitness(weights);
    sizeBitVector_ = (vertices_.size() / INT_SIZE) + 1;

    //clique set
    cliques_ = new CoinCliqueList(4096, 32768);
    clqWeightCap_ = 4096;
    clqWeight_ = std::vector<double>(clqWeightCap_);

    nC_ = 0;
    weightC_ = 0.0;
    C_ = std::vector<size_t>(maxDegree + 1);

    allIn_ = std::vector<size_t>(sizeBitVector_);

    nS_ = std::vector<size_t>(maxDegree + 2);
    nP_ = std::vector<size_t>(maxDegree + 2);

    S_ = std::vector<std::vector<size_t> >(maxDegree + 2, std::vector<size_t>(sizeBitVector_));
    L_ = std::vector<std::vector<size_t> >(maxDegree + 2, std::vector<size_t>(sizeBitVector_));
    P_ = std::vector<std::vector<size_t> >(maxDegree + 2, std::vector<size_t>(sizeBitVector_));

    mask_ = std::vector<size_t>(INT_SIZE);
    mask_[0] = 1;
    for (size_t h = 1; h < INT_SIZE; h++) {
        mask_[h] = mask_[h - 1] << 1U;
    }

    cgBitstring_ = std::vector<std::vector<size_t> >(vertices_.size(), std::vector<size_t>(sizeBitVector_));
    ccgBitstring_ = std::vector<std::vector<size_t> >(vertices_.size(), std::vector<size_t>(sizeBitVector_));

    for(size_t u = 0; u < vertices_.size(); u++) {
        BKVertex &vertexU = vertices_[u];

        allIn_[u / INT_SIZE] |= mask_[u % INT_SIZE];
        ccgBitstring_[u][u / INT_SIZE] |= mask_[u % INT_SIZE];

        for(size_t v = u + 1; v < vertices_.size(); v++) {
            BKVertex &vertexV = vertices_[v];

            if (cgraph_->conflicting(vertexU.idx, vertexV.idx)) {
                cgBitstring_[u][v / INT_SIZE] |= mask_[v % INT_SIZE];
                cgBitstring_[v][u / INT_SIZE] |= mask_[u % INT_SIZE];
            } else {
                ccgBitstring_[u][v / INT_SIZE] |= mask_[v % INT_SIZE];
                ccgBitstring_[v][u / INT_SIZE] |= mask_[u % INT_SIZE];
            }
        }
    }
}


CoinBronKerbosch::~CoinBronKerbosch() {
    if (!vertices_.empty()) {
        delete cliques_;
    }
}

double CoinBronKerbosch::weightP(const size_t depth, size_t &u) {
    double weight = 0.0;
    bool first = false;

    for (size_t t = 0; t < sizeBitVector_; t++) {
        size_t value = P_[depth][t];

        while (value) {
            const size_t node = ((size_t)log2(value & (~value + 1u))) + (INT_SIZE * t);
            weight += vertices_[node].weight;
            value &= value - 1;

            if (!first) {
                u = node;
                first = true;
            }
        }
    }

#ifdef DEBUGCG
    assert(first);
#endif

    return weight;
}

void CoinBronKerbosch::bronKerbosch(size_t depth) {
    //P and S are empty
    if( (nP_[depth] == 0) && (nS_[depth] == 0) && (nC_ > 0) && (weightC_ >= minWeight_) ) {
        const size_t nCliques = cliques_->nCliques();
        //checking memory
        if (nCliques + 1 > clqWeightCap_) {
            clqWeightCap_ *= 2;
            clqWeight_.resize(clqWeightCap_);
        }
        //maximal clique above a threshold found
        cliques_->addClique(nC_, C_.data());
        clqWeight_[nCliques] = weightC_;

#ifdef DEBUGCG
        CoinCliqueList::validateClique(cgraph_, C_.data(), nC_);
#endif

        return;
    }

    if(calls_ > maxCalls_) {
        completeSearch_ = false;
        return;
    }

    if (nP_[depth] == 0) {
        return;
    }

    size_t u; //node with the highest fitness
    const double wP = weightP(depth, u);

    if(weightC_ + wP >= minWeight_) {
        //L = P \ N(u)
        for(size_t i = 0; i < sizeBitVector_; i++) {
            L_[depth][i] = P_[depth][i] & ccgBitstring_[u][i];
        }

        //for each v in L
        for (size_t t = 0; t < sizeBitVector_; t++) {
            size_t value = L_[depth][t];

            while (value) {
                const size_t v = ((size_t)log2(value & (~value + 1u))) + (INT_SIZE * t);
                value &= value - 1;

                //new P and S
                nS_[depth + 1] = nP_[depth + 1] = 0;
                for (size_t i = 0; i < sizeBitVector_; i++) {
                    S_[depth + 1][i] = S_[depth][i] & cgBitstring_[v][i];
                    P_[depth + 1][i] = P_[depth][i] & cgBitstring_[v][i];

                    if (S_[depth + 1][i]) {
                        nS_[depth + 1]++;
                    }

                    if (P_[depth + 1][i]) {
                        nP_[depth + 1]++;
                    }
                }

                //new C
                C_[nC_++] = vertices_[v].idx;
                weightC_ += vertices_[v].weight;

                //recursive call
                calls_++;
                bronKerbosch(depth + 1);

                //restoring C
                nC_--;
                weightC_ -= vertices_[v].weight;

                //P = P \ {v}
                const size_t backup = allIn_[v / INT_SIZE];
                allIn_[v / INT_SIZE] = allIn_[v / INT_SIZE] & (~mask_[v % INT_SIZE]);
                assert(nP_[depth] > 0);
                nP_[depth] = 0;
                for (size_t i = 0; i < sizeBitVector_; i++) {
                    P_[depth][i] = P_[depth][i] & allIn_[i];

                    if (P_[depth][i]) {
                        nP_[depth]++;
                    }
                }
                allIn_[v / INT_SIZE] = backup;

                //S = S U {v}
                S_[depth][v / INT_SIZE] |= mask_[v % INT_SIZE];//adding v in S
                nS_[depth]--;
            }
        }
    }
}

void CoinBronKerbosch::findCliques() {
    if (vertices_.empty()) {
        return;
    }

    nP_[0] = vertices_.size();
    for(size_t i = 0; i < vertices_.size(); i++) {
        P_[0][i / INT_SIZE] |= mask_[i % INT_SIZE];
    }

    nC_ = calls_ = 0;
    weightC_ = 0.0;
    completeSearch_ = true;
    bronKerbosch(0);
}

size_t CoinBronKerbosch::nCliques() const {
    if (!cliques_) {
        return 0;
    }
    return cliques_->nCliques();
}

const size_t* CoinBronKerbosch::getClique(const size_t i) const {
#ifdef DEBUGCG
    assert(i < cliques_->nCliques());
#endif
    return cliques_->cliqueElements(i);
}

size_t CoinBronKerbosch::getCliqueSize(const size_t i) const {
#ifdef DEBUGCG
    assert(i < cliques_->nCliques());
#endif
    return cliques_->cliqueSize(i);
}

double CoinBronKerbosch::getCliqueWeight(const size_t i) const {
#ifdef DEBUGCG
    assert(i < cliques_->nCliques());
#endif
    return clqWeight_[i];
}

void CoinBronKerbosch::setMinWeight(double minWeight) {
    minWeight_ = minWeight;
}

void CoinBronKerbosch::setMaxCalls(size_t maxCalls) {
    maxCalls_ = maxCalls;
}

bool CoinBronKerbosch::completedSearch() const {
    return completeSearch_;
}

size_t CoinBronKerbosch::numCalls() const {
    return calls_;
}

void CoinBronKerbosch::computeFitness(const double *weights) {
    switch (pivotingStrategy_) {
        case PivotingStrategy::Off:
            //do nothing
            break;
        case PivotingStrategy::Random: {
            shuffle_vertices(vertices_.data(), vertices_.size());
            break;
        }
        case PivotingStrategy::Degree: {
            for (size_t u = 0; u < vertices_.size(); u++) {
                const size_t uIdx = vertices_[u].idx;
                vertices_[u].fitness = cgraph_->degree(uIdx);
            }
            std::sort(vertices_.begin(), vertices_.end(), compareNodes);
            break;
        }
        case PivotingStrategy::Weight: {
            for (size_t u = 0; u < vertices_.size(); u++) {
                const size_t uIdx = vertices_[u].idx;
                vertices_[u].fitness = weights[uIdx];
            }
            std::sort(vertices_.begin(), vertices_.end(), compareNodes);
            break;
        }
        case PivotingStrategy::ModifiedDegree: {
            for (size_t u = 0; u < vertices_.size(); u++) {
                const size_t uIdx = vertices_[u].idx;
                vertices_[u].fitness = cgraph_->modifiedDegree(uIdx);
            }
            std::sort(vertices_.begin(), vertices_.end(), compareNodes);
            break;
        }
        case PivotingStrategy::ModifiedWeight: {
            std::vector<size_t> neighs = std::vector<size_t>(cgraph_->size());
            std::vector<char> iv = std::vector<char>(cgraph_->size());
            for (size_t u = 0; u < vertices_.size(); u++) {
                const size_t uIdx = vertices_[u].idx;
                const std::pair<size_t, const size_t*> rescg = cgraph_->conflictingNodes(uIdx, neighs.data(), iv.data());
                vertices_[u].fitness = weights[uIdx];
                for (size_t v = 0; v < rescg.first; v++) {
                    const size_t vIdx = rescg.second[v];
                    vertices_[u].fitness += weights[vIdx];
                }
            }
            std::sort(vertices_.begin(), vertices_.end(), compareNodes);
            break;
        }
        case PivotingStrategy::ModifiedDegreeWeight: {
            std::vector<size_t> neighs = std::vector<size_t>(cgraph_->size());
            std::vector<char> iv = std::vector<char>(cgraph_->size());
            for (size_t u = 0; u < vertices_.size(); u++) {
                const size_t uIdx = vertices_[u].idx;
                const std::pair<size_t, const size_t*> rescg = cgraph_->conflictingNodes(uIdx, neighs.data(), iv.data());
                vertices_[u].fitness = weights[uIdx] + cgraph_->modifiedDegree(uIdx);
                for (size_t v = 0; v < rescg.first; v++) {
                    const size_t vIdx = rescg.second[v];
                    vertices_[u].fitness += (weights[vIdx] + cgraph_->modifiedDegree(vIdx));
                }
            }
            std::sort(vertices_.begin(), vertices_.end(), compareNodes);
            break;
        }
        default:
            fprintf(stderr, "Invalid option %u for pivoting strategy!\n", pivotingStrategy_);
            abort();
    }
}

static void shuffle_vertices (BKVertex *vertices, size_t n) {
    if (n <= 1) {
        return;
    }

    BKVertex tmp;

    srand (time(NULL));

    for (size_t i = n - 1; i > 0; i--) {
        size_t j = rand() % (i + 1);
        tmp = vertices[i];
        vertices[i] = vertices[j];
        vertices[j] = tmp;
    }
}
