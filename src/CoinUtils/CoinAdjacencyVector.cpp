/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class to store a (growable) list of neighbors for each node
 * Initially implemented to be used in the Conflict Graph
 *
 * @file CoinAdjacencyVector.cpp
 * @brief Vector of growable vectors
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <cassert>
#include <algorithm>

#include "CoinAdjacencyVector.hpp"

#define KEEP_SORTED

CoinAdjacencyVector::CoinAdjacencyVector( size_t _nRows, size_t _iniRowSize )
  : rows_(std::vector<std::vector<size_t> >(_nRows))
  , notUpdated_(std::vector<size_t>(_nRows, 0))
{
    for (size_t i=0; i < rows_.size(); ++i)
        rows_[i].reserve(_iniRowSize);
}

CoinAdjacencyVector::~CoinAdjacencyVector()
{
}

const size_t *CoinAdjacencyVector::getRow ( size_t idxRow ) const
{
  assert(idxRow<this->rows_.size());

  return rows_[idxRow].data();
}

bool CoinAdjacencyVector::isNeighbor(size_t idxNode, size_t idxNeigh) const {
#ifdef KEEP_SORTED
    return std::binary_search(rows_[idxNode].begin(), rows_[idxNode].end(), idxNeigh);
#else
    for (std::vector<size_t>::const_iterator i=rows_[idxNode].begin(); i != rows_[idxNode].end(); ++i)
        if (*i == idxNeigh)
            return true;
#endif

    return false;
}

void CoinAdjacencyVector::addNeighbor( size_t idxNode, size_t idxNeigh, bool addReverse ) {
  checkCapNode(idxNode);
#ifdef KEEP_SORTED
  bool res = tryAddElementSortedVector(rows_[idxNode], idxNeigh);
  if (res) {
    if (addReverse)
      addNeighbor( idxNeigh, idxNode, false );
  }
#else
    rows_[idxNode].push_back(idxNeigh);
    if (addReverse)
      this->addNeighbor(idxNeigh, idxNode, false);
#endif
}


size_t CoinAdjacencyVector::rowSize( size_t idxRow ) const
{
  assert(idxRow<this->rows_.size());

  return rows_[idxRow].size();
}

void CoinAdjacencyVector::checkCapNode( const size_t idxNode, const size_t newEl )
{
    assert( idxNode < rows_.size() );

    size_t currCap = rows_[idxNode].capacity();
    size_t currSize = rows_[idxNode].size();

    // no need to resize
    if ( currSize + newEl <= currCap )
      return;

    // for resizing
    const size_t newIdxNodeCap = std::max(currCap*2, currSize+newEl);
    rows_[idxNode].reserve(newIdxNodeCap);
}

void CoinAdjacencyVector::fastAddNeighbor( size_t idxNode, size_t idxNeigh )
{
  //printf("adding to %zu %zu currCap: %zu currSize: %zu\n", idxNode, idxNeigh, rows_[idxNode].capacity(), rows_[idxNode].size()); fflush(stdout);

  checkCapNode(idxNode);

  rows_[idxNode].push_back(idxNeigh);
  notUpdated_[idxNode]++;
}

void CoinAdjacencyVector::sort()
{
  for ( size_t i=0 ; i<rows_.size() ; ++i )
    std::sort(rows_[i].begin(), rows_[i].end());
}

bool CoinAdjacencyVector::tryAddElementSortedVector(std::vector<size_t> &el, size_t newEl)
{
  if (el.size()) {
    std::vector<size_t>::iterator pos = std::lower_bound(el.begin(), el.end(), newEl);
    if (pos != el.end() && *pos == newEl)
      return false;
    el.insert(pos, newEl);
  } else {
    el.push_back(newEl);
  }
  return true;
}

size_t CoinAdjacencyVector::totalElements() const
{
  size_t res = 0;

  for ( size_t i=0 ; i<rows_.size() ; ++i )
    res += rows_[i].size();

  return res;
}

void CoinAdjacencyVector::flush() {
  for ( size_t i=0 ; i<this->rows_.size() ; ++i ) {
    if (notUpdated_[i]) {
      std::sort(rows_[i].begin(), rows_[i].end());
      std::vector<size_t>::iterator newEnd = std::unique(rows_[i].begin(), rows_[i].end());
      rows_[i].erase(newEnd, rows_[i].end());
      notUpdated_[i] = 0;
    }
  }
}

void CoinAdjacencyVector::addNeighborsBuffer( size_t idxNode, size_t n, const size_t elements[] ) {
  checkCapNode(idxNode, n);
  for ( size_t i=0 ; (i<n) ; ++i )
    if (elements[i] != idxNode) {
      rows_[idxNode].push_back(elements[i]);
      notUpdated_[idxNode]++;
    }
}

void CoinAdjacencyVector::sort(size_t idxRow) {
  std::sort(rows_[idxRow].begin(), rows_[idxRow].end());
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
