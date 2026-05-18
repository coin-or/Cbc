/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class that implements Dijkstra's algorithm for
 * finding the shortest paths between nodes in a graph.
 * Specialized for sparse graphs.
 *
 * @file CoinShortestPath.cpp
 * @brief Shortest path algorithm
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include "CoinShortestPath.hpp"
#include "CoinNodeHeap.hpp"

#define SPATH_EPS 1e-6
#define SPATH_INFTY_NODE std::numeric_limits<size_t>::max()
#define SPATH_INFTY_DIST std::numeric_limits<double>::max()

CoinShortestPath::CoinShortestPath(size_t nodes, size_t arcs, const size_t *arcStart, const size_t *toNode, const double *dist) {
    nodes_ = nodes;
    nh_ = new CoinNodeHeap(nodes_);
    previous_ = std::vector<size_t>(nodes_);
    dist_ = std::vector<double>(nodes_);
    path_ = std::vector<size_t>(nodes_);
    neighs_ = std::vector<std::vector<std::pair<size_t, double> > >(nodes_);

    size_t idx = 0;
    for (size_t n = 0; n < nodes_; n++) {
        neighs_[n] = std::vector<std::pair<size_t, double> >(arcStart[n+1] - arcStart[n]);
        for (size_t i = 0; i < neighs_[n].size(); ++i) {
            neighs_[n][i].first = toNode[idx];
            neighs_[n][i].second = dist[idx];
#ifdef DEBUGCG
            assert(neighs_[n][i].first < nodes);
#endif
            idx++;
        }
    }
}

CoinShortestPath::~CoinShortestPath() {
    delete nh_;
}

void CoinShortestPath::find(const size_t origin) {
    assert(origin < this->nodes_);

    nh_->reset();

    for (size_t i = 0; i < nodes_; i++) {
        previous_[i] = SPATH_INFTY_NODE;
        dist_[i] = SPATH_INFTY_DIST;
    }

    dist_[origin] = 0.0;
    nh_->update(origin, 0.0);

    size_t topNode;
    double topCost;
    while (!nh_->isEmpty()) {
        topCost = nh_->removeFirst(&topNode);
#ifdef DEBUGCG
        assert(topCost + SPATH_EPS <= SPATH_INFTY_DIST);
#endif
        // updating neighbors distances by iterating in all neighbors
        for (std::vector<std::pair<size_t, double> >::iterator n = neighs_[topNode].begin(); n != neighs_[topNode].end(); ++n) {
            const size_t toNode = n->first;
            const double dist = n->second;
            const double newDist = topCost + dist;

            assert( dist >= 0 );
            assert( newDist >= 0 );

            if (dist_[toNode] >= newDist + SPATH_EPS) {
                previous_[toNode] = topNode;
                dist_[toNode] = newDist;
                nh_->update(toNode, newDist);
            } // updating heap if necessary
        } // going through node neighbors
    } // going through all nodes in priority queue
}

void CoinShortestPath::find(const size_t origin, const size_t destination) {
    assert(origin < this->nodes_);
    assert(destination < this->nodes_);
    nh_->reset();

    for (size_t i = 0; i < nodes_; i++) {
        previous_[i] = SPATH_INFTY_NODE;
        dist_[i] = SPATH_INFTY_DIST;
    }

    dist_[origin] = 0.0;
    nh_->update(origin, 0.0);

    size_t topNode;
    double topCost;
    while (!nh_->isEmpty()) {
        topCost = nh_->removeFirst(&topNode);
#ifdef DEBUGCG
        assert(topCost + SPATH_EPS <= SPATH_INFTY_DIST);
#endif
        if(topNode == destination) {
            break;
        }
        // updating neighbors distances by iterating in all neighbors
        for (std::vector<std::pair<size_t, double> >::iterator n = neighs_[topNode].begin(); n != neighs_[topNode].end(); ++n) {
            const size_t toNode = n->first;
            const double dist = n->second;
            const double newDist = topCost + dist;

            if (dist_[toNode] >= newDist + SPATH_EPS) {
                previous_[toNode] = topNode;
                dist_[toNode] = newDist;
                nh_->update(toNode, newDist);
            } // updating heap if necessary
        } // going through node neighbors
    } // going through all nodes in priority queue
}

size_t CoinShortestPath::path(size_t toNode, size_t *indexes) {
    assert(toNode < this->nodes_);

    size_t n = 0, currNode = toNode;

    path_[n++] = currNode;

    while((currNode = previous_[currNode]) != SPATH_INFTY_NODE) {
        path_[n++] = currNode;
    }

    for (size_t i = 0; i < n; i++) {
        indexes[i] = path_[n-i-1];
    }

    return n;
}

double CoinShortestPath::distance(size_t node) const {
#ifdef DEBUGCG
    assert(node < nodes_);
#endif
    return dist_[node];
}

size_t CoinShortestPath::previous(size_t node) const {
#ifdef DEBUGCG
    assert(node < nodes_);
#endif
    return previous_[node];
}

