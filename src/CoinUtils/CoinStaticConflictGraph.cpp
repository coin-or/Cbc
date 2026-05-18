/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * CoinConflictGraph implementation which supports fast queries
 * but doesn't supports modifications.
 *
 * @file CoinStaticConflictGraph.cpp
 * @brief static CoinConflictGraph implementation with fast queries
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <algorithm>
#include <cstring>
#include <limits>
#include <cstdlib>
#include <cstdio>

#include "CoinStaticConflictGraph.hpp"
#include "CoinDynamicConflictGraph.hpp"
#include "CoinCliqueList.hpp"

CoinStaticConflictGraph::CoinStaticConflictGraph ( const CoinConflictGraph *cgraph )
{
  iniCoinStaticConflictGraph(cgraph);
}

CoinStaticConflictGraph::CoinStaticConflictGraph (
  const int numCols,
  const char *colType,
  const double *colLB,
  const double *colUB,
  const CoinPackedMatrix *matrixByRow,
  const char *sense,
  const double *rowRHS,
  const double *rowRange,
  const double primalTolerance,
  const double infinity,
  const std::vector<std::string> &colNames,
  const std::vector<std::string> &rowNames)
{
    CoinDynamicConflictGraph *cgraph = new CoinDynamicConflictGraph(numCols, colType, colLB, colUB, matrixByRow, sense, rowRHS, rowRange, primalTolerance, infinity, colNames, rowNames);

    iniCoinConflictGraph(cgraph);
    nDirectConflicts_ = cgraph->nTotalDirectConflicts();
    totalCliqueElements_ = cgraph->nTotalCliqueElements();

    degree_ = std::vector<size_t>(size_, 0);
    modifiedDegree_ = std::vector<size_t>(size_, 0);
    conflicts_ = std::vector<std::vector<size_t>>(size_);
    nodeCliques_ = std::vector<std::vector<size_t>>(size_);
    cliques_ = std::vector<std::vector<size_t>>(cgraph->nCliques());
    infeasibleImplications_ = cgraph->infeasibleImplications();

    // move direct conflicts from dynamic graph (avoids copy)
    for (size_t i = 0; i < size(); ++i)
      conflicts_[i] = cgraph->moveDirectConflicts(i);

    // move cliques from dynamic graph
    for (size_t ic = 0; ic < cgraph->nCliques(); ++ic)
      cliques_[ic] = cgraph->moveClique(ic);

    // filling node cliques
    for (size_t ic = 0; ic < cliques_.size(); ++ic) {
      const size_t clqSize = cliques_[ic].size();
      for (size_t j = 0; j < clqSize; ++j)
        nodeCliques_[cliques_[ic][j]].push_back(ic);
    }

    for (size_t i = 0; i < size_; ++i) {
      this->setDegree(i, cgraph->degree(i));
      this->setModifiedDegree(i, cgraph->modifiedDegree(i));
    }

    newBounds_ = cgraph->updatedBounds();

#ifdef CGRAPH_STATS
    memcpy(rowTypeStats_, cgraph->rowTypeStats(), sizeof(rowTypeStats_));
    rowProfileStats_ = cgraph->rowProfileStats();
#endif

    delete cgraph;
}

bool CoinStaticConflictGraph::nodeInClique( size_t idxClique, size_t node ) const
{
  return std::binary_search(cliques_[idxClique].begin(), cliques_[idxClique].end(), node);
}

CoinStaticConflictGraph *CoinStaticConflictGraph::clone() const
{
  return new CoinStaticConflictGraph ( this );
}

size_t CoinStaticConflictGraph::nDirectConflicts ( size_t idxNode ) const
{
  return this->conflicts_[idxNode].size();
}

const size_t * CoinStaticConflictGraph::directConflicts ( size_t idxNode ) const
{
  return this->conflicts_[idxNode].data();
}

size_t CoinStaticConflictGraph::nCliques() const
{
  return this->cliques_.size();
}

const size_t * CoinStaticConflictGraph::cliqueElements ( size_t idxClique ) const
{
  return this->cliques_[idxClique].data();
}

size_t CoinStaticConflictGraph::cliqueSize( size_t idxClique ) const {
  return this->cliques_[idxClique].size();
}

const size_t * CoinStaticConflictGraph::nodeCliques ( size_t idxNode ) const
{
  return nodeCliques_[idxNode].data();
}

size_t CoinStaticConflictGraph::nNodeCliques ( size_t idxNode ) const
{
  return this->nodeCliques_[idxNode].size();
}

void CoinStaticConflictGraph::setDegree(size_t idxNode, size_t deg)
{
  this->degree_[idxNode] = deg;
}

void CoinStaticConflictGraph::setModifiedDegree(size_t idxNode, size_t mdegree)
{
    this->modifiedDegree_[idxNode] = mdegree;
}

size_t CoinStaticConflictGraph::degree(const size_t node) const
{
  return degree_[node];
}

size_t CoinStaticConflictGraph::modifiedDegree(const size_t node) const
{
    return modifiedDegree_[node];
}

CoinStaticConflictGraph::CoinStaticConflictGraph( const CoinConflictGraph *cgraph, const size_t n, const size_t elements[] )
{
  iniCoinConflictGraph( n );

#define REMOVED std::numeric_limits< size_t >::max()
  nDirectConflicts_ = totalCliqueElements_ = 0;

  std::vector< size_t > newIdx( cgraph->size(), REMOVED );
  for ( size_t i=0 ; (i<n) ; ++i )
    newIdx[elements[i]] = i;

  std::vector<char> iv(size_);
  std::vector< bool > ivNeighs;

  // large and small cliques set
  CoinCliqueList smallClqs( 4096, 32768 );
  CoinCliqueList largeClqs( 4096, 32768 );

  std::vector<size_t> clqEls(size_);

  // separating new cliques (removing variables) into small and large
  for ( size_t ic = 0 ; (ic<cgraph->nCliques()) ; ++ic ) {
    size_t nEl = 0;
    for ( size_t j=0 ; (j<cgraph->cliqueSize(ic)) ; ++j ) {
      size_t idxNode = newIdx[ cgraph->cliqueElements(ic)[j] ];
      if ( idxNode == REMOVED )
        continue;

      clqEls[nEl++] = idxNode;
    }

    if ( nEl >= CoinConflictGraph::minClqRow_ ) {
      largeClqs.addClique( nEl, clqEls.data() );
    } else {
      smallClqs.addClique( nEl, clqEls.data() );
    }
  }

//  printf("In induced subgraph there are still %zu large cliques and %zu cliques will now be stored as pairwise conflicts.\n",
//    largeClqs.nCliques(), smallClqs.nCliques() ); fflush( stdout );

  // checking in small cliques new direct neighbors of each node
  CoinAdjacencyVector newNeigh( size_, 16 );

  smallClqs.computeNodeOccurrences( this->size() );
  largeClqs.computeNodeOccurrences( this->size() );

  // computing new direct conflicts
  for ( size_t i=0 ; (i<smallClqs.nDifferentNodes()) ; ++i ) {
    size_t idxNode = smallClqs.differentNodes()[i];

    iv[idxNode] = 1;

    size_t idxOrigNode = elements[idxNode];

    // marking known direct conflicts
    for ( size_t j=0 ; (j<cgraph->nDirectConflicts(idxOrigNode)) ; ++j )
      if ( newIdx[cgraph->directConflicts(idxOrigNode)[j]] != REMOVED )
        iv[newIdx[cgraph->directConflicts(idxOrigNode)[j]]] = 1;

    // marking those that appear in the large cliques
    for ( size_t j=0 ; j<largeClqs.nNodeOccurrences(idxNode) ; ++j ) {
      size_t idxClq = largeClqs.nodeOccurrences(idxNode)[j];

      // all elements of this large clique
      for ( size_t j=0 ; (j<largeClqs.cliqueSize(idxClq)) ; ++j )
        iv[largeClqs.cliqueElements(idxClq)[j]] = 1;
    }

    // checking with neighbors from small cliques are not
    // yet in direct conflicts or in the remaining large cliques
    for ( size_t j=0 ; (j<smallClqs.nNodeOccurrences(idxNode)) ; ++j ) {
      size_t idxClq = smallClqs.nodeOccurrences(idxNode)[j];
      for ( size_t k=0 ; (k<smallClqs.cliqueSize(idxClq)) ; ++k ) {
        if (!iv[smallClqs.cliqueElements(idxClq)[k]]) {
          iv[smallClqs.cliqueElements(idxClq)[k]] = 1;
          newNeigh.fastAddNeighbor( idxNode, smallClqs.cliqueElements(idxClq)[k] );
        }
      }
    }

    newNeigh.flush();

    // marking know direct conflicts
    for ( size_t j=0 ; (j<cgraph->nDirectConflicts(idxOrigNode)) ; ++j )
      if ( newIdx[cgraph->directConflicts(idxOrigNode)[j]] != REMOVED )
        iv[newIdx[cgraph->directConflicts(idxOrigNode)[j]]] = 0;

    // marking those that appear in the large cliques
    for ( size_t j=0 ; j<largeClqs.nNodeOccurrences(idxNode) ; ++j ) {
      size_t idxClq = largeClqs.nodeOccurrences(idxNode)[j];

      // all elements of this large clique
      for ( size_t j=0 ; (j<largeClqs.cliqueSize(idxClq)) ; ++j )
        iv[largeClqs.cliqueElements(idxClq)[j]] = 0;
    }

    // unchecking new direct conflicts
    for ( size_t j=0 ; (j<newNeigh.rowSize(idxNode)) ; ++j )
      iv[newNeigh.getRow(idxNode)[j]] = 0;

    iv[idxNode] = 0;
  }

  // computing new number of direct conflicts per node
  size_t prevTotalDC = 0;
  std::vector<size_t> prevDC(size_);

  for ( size_t i=0 ; (i<n) ; ++i ) {
    size_t idxOrig = elements[i];
    prevDC[i] = 0;

    for ( size_t j=0 ; ( j < cgraph->nDirectConflicts(idxOrig) ) ; ++j ) {
      size_t ni = newIdx[ cgraph->directConflicts(idxOrig)[j] ] ;
      if ( ni == REMOVED )
        continue;
      prevDC[i]++;
    }

    prevTotalDC += prevDC[i];
  }

  nDirectConflicts_ = prevTotalDC + newNeigh.totalElements();
  totalCliqueElements_ = largeClqs.totalElements();
  degree_ = std::vector<size_t>(size_);
  modifiedDegree_ = std::vector<size_t>(size_);
  conflicts_ = std::vector<std::vector<size_t> >(size_);
  nodeCliques_ = std::vector<std::vector<size_t> >(size_);
  cliques_ = std::vector<std::vector<size_t> >(largeClqs.nCliques());

  // filling cliques
  for ( size_t i=0 ; (i<largeClqs.nCliques()) ; ++i ) {
    cliques_[i] = std::vector<size_t>(largeClqs.cliqueElements(i), largeClqs.cliqueElements(i) + largeClqs.cliqueSize(i));
  }


  // copying remaining direct conflicts
  // adding new conflicts when they exist
  for ( size_t i=0 ; (i<n) ; ++i ) {
    size_t idxOrig = elements[i];
    std::vector<size_t> conf;

    for ( size_t j=0 ; ( j < cgraph->nDirectConflicts(idxOrig) ) ; ++j ) {
      size_t ni = newIdx[ cgraph->directConflicts(idxOrig)[j] ] ;
      if ( ni == REMOVED )
        continue;
      conf.push_back(ni);
    }

    conflicts_[i] = conf;

    // new pairwise conflicts from new smallCliques
    if (newNeigh.rowSize(i)) {
      conflicts_[i].insert(conflicts_[i].end(), newNeigh.getRow(i), newNeigh.getRow(i) + newNeigh.rowSize(i));
      std::sort(conflicts_[i].begin(), conflicts_[i].end());
    }
  } // all nodes

  // filling node cliques
  for ( size_t i=0 ; i<size_ ; ++i )
    if (largeClqs.nNodeOccurrences(i))
      nodeCliques_[i] = std::vector<size_t>(largeClqs.nodeOccurrences(i), largeClqs.nodeOccurrences(i) + largeClqs.nNodeOccurrences(i));

  this->recomputeDegree();
#undef REMOVED
}

size_t CoinStaticConflictGraph::nTotalDirectConflicts() const {
  return this->nDirectConflicts_;
}

size_t CoinStaticConflictGraph::nTotalCliqueElements() const {
  return this->totalCliqueElements_;
}

CoinStaticConflictGraph::~CoinStaticConflictGraph()
{
}

const std::vector< std::pair< size_t, std::pair< double, double > > > & CoinStaticConflictGraph::updatedBounds() const
{
    return newBounds_;
}

void CoinStaticConflictGraph::iniCoinStaticConflictGraph(const CoinConflictGraph *cgraph) {
    iniCoinConflictGraph(cgraph);
    nDirectConflicts_  = cgraph->nTotalDirectConflicts();
    totalCliqueElements_ = cgraph->nTotalCliqueElements();

    degree_ = std::vector<size_t>(size_, 0);
    modifiedDegree_ = std::vector<size_t>(size_, 0);
    conflicts_ = std::vector<std::vector<size_t> >(size_);
    nodeCliques_ = std::vector<std::vector<size_t> >(size_);
    cliques_ = std::vector<std::vector<size_t> >(cgraph->nCliques());
    infeasibleImplications_ = cgraph->infeasibleImplications();

    // copying direct conflicts
    for ( size_t i=0 ; (i<size()) ; ++i ) {
        const size_t sizeConf = cgraph->nDirectConflicts(i);
        const size_t *conf = cgraph->directConflicts(i);
        conflicts_[i] = std::vector<size_t>(conf, conf + sizeConf);
    } // all nodes

    // copying cliques
    for ( size_t ic=0 ; ( ic<(size_t)cgraph->nCliques() ) ; ++ic )
    {
        const size_t *clique = cgraph->cliqueElements(ic);
        const size_t *cliqueEnd = clique + cgraph->cliqueSize(ic);
        // copying clique contents
        cliques_[ic] = std::vector<size_t>(clique, cliqueEnd);
    }

    // filling node cliques
    for ( size_t ic=0 ; ( ic < cgraph->nCliques() ) ; ++ic )
    {
        const size_t *clq = cliqueElements(ic);
        const size_t clqSize = cliques_[ic].size();
        for ( size_t iclqe=0 ; (iclqe<clqSize) ; ++iclqe )
        {
            size_t el = clq[iclqe];
            nodeCliques_[el].push_back(ic);
        }
    }

    for ( size_t i=0 ; (i<size_) ; ++i ) {
        this->setDegree(i, cgraph->degree(i));
        this->setModifiedDegree(i, cgraph->modifiedDegree(i));
    }
}
