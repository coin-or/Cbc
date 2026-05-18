/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * A class to store a sequential list of cliques.
 * Optionally computes occurrence of nodes in cliques.
 *
 * @file CoinCliqueList.hpp
 * @brief List of cliques
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef COINCLIQUELIST_H
#define COINCLIQUELIST_H

#include "CoinUtilsConfig.h"
#include <cstddef>
#include <vector>

#ifdef DEBUGCG
    class CoinConflictGraph;
#endif

/**
 * A class to store a sequential list of cliques.
 * Optionally computes occurrence of nodes in cliques.
 */
class COINUTILSLIB_EXPORT CoinCliqueList
{
public:
  /**
   * Default constructor
   *
   * @param _iniClqCap initial capacity to store cliques
   * @param _iniClqElCap initial capacity of the elements
   * of the cliques
   **/
  CoinCliqueList( size_t _iniClqCap, size_t _iniClqElCap );

  /**
   * Destructor
   **/
  virtual ~CoinCliqueList();

  /**
   * Insert a clique in the list.
   **/
  void addClique( size_t size, const size_t els[] );

  /**
   * Return the number of cliques stored in the list.
   **/
  size_t nCliques() const;

  /**
   * Return the size of a clique.
   *
   * @param idxClq index of a clique.
   **/
  size_t cliqueSize( size_t idxClq ) const;

  /**
   * Return the elements of a clique.
   *
   * @param idxClq index of a clique.
   **/
  const size_t *cliqueElements( size_t idxClq ) const;

  /**
   * Move the contents of a clique out (transfers ownership, leaves it empty).
   **/
  std::vector<size_t> moveClique(size_t idxClq) { return std::move(clqEls_[idxClq]); }

  /**
   * Return the total number of elements
   * considering all cliques.
   **/
  size_t totalElements() const;

  /**
   * Fill the data structures indicating
   * in which clique each node appears.
   **/
  void computeNodeOccurrences( size_t nNodes );

  /**
   * Return the number of cliques in which
   * a node appears.
   *
   * @param idxNode index of the node
   **/
  size_t nNodeOccurrences( size_t idxNode ) const;

  /**
   * Return the indexes of the cliques in which
   * a node appears.
   *
   * @param idxNode index of the node
   **/
  const size_t *nodeOccurrences( size_t idxNode) const;

  /**
   * Return the number of different nodes
   * that are stored in the clique list.
   **/
  size_t nDifferentNodes() const;

  /**
   * Return the nodes (without duplicates)
   * that are stored in the clique list.
   **/
  const size_t *differentNodes() const;

#ifdef DEBUGCG
    /**
     * Validate a clique according to a conflict graph.
     * returns false if bad in some way
     *
     * @param cgraph conflict graph
     * @param idxs indexes of the clique
     * @param size size of the clique
     **/
    static bool validateClique(const CoinConflictGraph *cgraph, const size_t *idxs, const size_t size);
#endif

private:
  /**
   * Capacity to store cliques in the list.
   **/
  size_t cliquesCap_;

  /**
   * Total number of elements considering all cliques.
   **/
  size_t nCliqueElements_;

  /**
   * Pointer with the elements of the cliques
   **/
  std::vector<std::vector<size_t> > clqEls_;

  /**
   * Store in which clique each node appears.
   * Only filled if computeNodeOccurrences is called.
   **/
  std::vector<std::vector<size_t> > nodeOccur_;

  /**
   * Number of different nodes that are stored
   * in the clique list.
   **/
  size_t nDifferent_;

  /**
   * Nodes (without duplicates) that are stored
   * in the clique list.
   **/
  std::vector<size_t> diffNodes_;
};

#endif // COINCLIQUELIST_H
