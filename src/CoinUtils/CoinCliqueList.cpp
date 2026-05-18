/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * A class to store a sequential list of cliques.
 * Optionally computes occurrence of nodes in cliques.
 *
 * @file CoinCliqueList.cpp
 * @brief List of cliques
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cassert>
#include "CoinCliqueList.hpp"

#ifdef DEBUGCG
    #include "CoinConflictGraph.hpp"
#endif

CoinCliqueList::CoinCliqueList( size_t _iniClqCap, size_t _iniClqElCap )
  : cliquesCap_( _iniClqCap )
  , nCliqueElements_( 0 )
  , nDifferent_ ( 0 )
{
}

void CoinCliqueList::addClique(size_t size, const size_t els[]) {
  if (clqEls_.size()+1 > cliquesCap_) {
    cliquesCap_ *= 2;
    this->clqEls_.reserve(cliquesCap_);
  }

  clqEls_.push_back(std::vector<size_t>(els, els + size));

  nCliqueElements_ += size;
}

size_t CoinCliqueList::cliqueSize( size_t idxClq ) const {
  return this->clqEls_[idxClq].size();
}

const size_t *CoinCliqueList::cliqueElements( size_t idxClq ) const {
  return this->clqEls_[idxClq].data();
}

CoinCliqueList::~CoinCliqueList()
{
}

void CoinCliqueList::computeNodeOccurrences(size_t nNodes)
{
  nodeOccur_ = std::vector<std::vector<size_t> >(nNodes);

  // couting number of occurrences for each node
  std::vector<size_t> noc(nNodes, 0);
  for (size_t i=0; i < clqEls_.size(); ++i) {
    for (size_t j=0; j < clqEls_[i].size(); ++j)
        noc[clqEls_[i][j]]++;
  }

  for (size_t in=0; in<nNodes; ++in)
    nodeOccur_[in] = std::vector<size_t>(noc[in]);

  std::fill(noc.begin(), noc.end(), 0);

  nDifferent_ = 0;
  for ( size_t ic=0 ; ic<nCliques() ; ++ic ) {
    for ( size_t j=0 ; (j<cliqueSize(ic)) ; ++j ) {
      size_t node = cliqueElements(ic)[j];
      if (!noc[node])
        nDifferent_++;
      nodeOccur_[node][noc[node]] = ic;
      ++noc[node];
    }
  }

  std::fill(noc.begin(), noc.end(), 0);

  diffNodes_ = std::vector<size_t>(nDifferent_);
  nDifferent_ = 0;
  for (size_t i=0; i < clqEls_.size(); ++i) {
    for (size_t j=0; j < clqEls_[i].size(); ++j) {
      size_t node = clqEls_[i][j];
      if (!noc[node])
        diffNodes_[nDifferent_++] = node;
      noc[node]++;
    }
  }
}

size_t CoinCliqueList::nCliques() const
{
  return clqEls_.size();
}

const size_t * CoinCliqueList::differentNodes() const
{
  return this->diffNodes_.data();
}

size_t CoinCliqueList::nDifferentNodes() const
{
  return this->nDifferent_;
}

const size_t * CoinCliqueList::nodeOccurrences(size_t idxNode) const
{
  return this->nodeOccur_[idxNode].data();
}

size_t CoinCliqueList::totalElements() const {
  return this->nCliqueElements_;
}

size_t CoinCliqueList::nNodeOccurrences(size_t idxNode) const
{
  return this->nodeOccur_[idxNode].size();
}

#ifdef DEBUGCG
bool CoinCliqueList::validateClique(const CoinConflictGraph *cgraph, const size_t *idxs, const size_t size) {
    if (size == 0) {
        fprintf(stderr, "Empty clique!\n");
        abort();
    }

    for (size_t j = 0; j < size; j++) {
        const size_t vidx = idxs[j];
        assert(vidx < cgraph->size());
    }

    for (size_t i = 0; i < size - 1; i++) {
        for (size_t j = i + 1; j < size; j++) {
            if ((!cgraph->conflicting(idxs[i], idxs[j])) || (idxs[i] == idxs[j])) {
                fprintf(stderr, "ERROR: Nodes %ld and %ld are not in conflict.\n", idxs[i], idxs[j]);
                return false;// may be infeasible abort();
            }
        }
    }
    return true;
}
#endif
