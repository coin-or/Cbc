/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Monotone heap
 * Updates MUST always decrease costs
 *
 * @file CoinNodeHeap.cpp
 * @brief Monotone heap
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include "CoinNodeHeap.hpp"
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <limits>

#define NODEHEAP_EPS 1e-6
#define NODEHEAP_INFTY std::numeric_limits<double>::max()

// position of root node in vector
#define rootPos(node) ((node > 0) ? ((((node + 1) / 2) - 1)) : std::numeric_limits<size_t>::max())

// position of the first child node in vector
#define childPos(node) ((node * 2) + 1)

CoinNodeHeap::CoinNodeHeap(size_t numNodes) {
#ifdef DEBUGCG
    assert(numNodes > 0);
#endif
    numNodes_ = numNodes;
    pq_ = std::vector<std::pair<size_t, double> >(numNodes);
    pos_ = std::vector<size_t>(numNodes);
    reset();
}

CoinNodeHeap::~CoinNodeHeap() {}

void CoinNodeHeap::reset() {
    for (size_t i = 0; i < numNodes_; i++) {
        pq_[i].first = i;
        pq_[i].second = NODEHEAP_INFTY;
        pos_[i] = i;
    }
}

void CoinNodeHeap::update(size_t node, double cost) {
    const size_t pos = pos_[node];
    size_t root, child = pos;

    assert(cost + NODEHEAP_EPS <= pq_[pos].second);
    pq_[pos].second = cost;

    while ((root = rootPos(child)) != std::numeric_limits<size_t>::max()) {
        if (pq_[root].second >= pq_[child].second + NODEHEAP_EPS) {
            std::swap(pq_[child], pq_[root]);
            pos_[pq_[root].first] = root;
            pos_[pq_[child].first] = child;
            child = root;
        } else {
            return;
        }
    }
}

double CoinNodeHeap::removeFirst(size_t *node) {
    const size_t posLastNode = numNodes_ - 1;
    double cost = pq_[0].second;

    (*node) = pq_[0].first;
    pq_[0] = pq_[posLastNode];
    pq_[posLastNode].first = (*node);
    pq_[posLastNode].second = NODEHEAP_INFTY;
    pos_[pq_[0].first] = 0;
    pos_[(*node)] = posLastNode;

    size_t root = 0;
    size_t child;
    while ((child = childPos(root)) < numNodes_) {
        // child with the smallest cost
        if ((child + 1 < numNodes_) && (pq_[child].second >= pq_[child + 1].second + NODEHEAP_EPS)) {
            child++;
        }

        if (pq_[root].second >= pq_[child].second + NODEHEAP_EPS) {
            std::swap(pq_[root], pq_[child]);
            pos_[pq_[root].first] = root;
            pos_[pq_[child].first] = child;
            root = child;
        } else {
            break;
        }
    }

    return cost;
}

bool CoinNodeHeap::isEmpty() const {
    return (pq_[0].second >= NODEHEAP_INFTY);
}

