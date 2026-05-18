/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Monotone heap
 * Updates MUST always decrease costs
 *
 * @file CoinNodeHeap.hpp
 * @brief Monotone heap
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef COINNODEHEAP_HPP
#define COINNODEHEAP_HPP

#include "CoinUtilsConfig.h"
#include <cstddef>
#include <utility>
#include <vector>

/**
 * Monotone heap.
 * Updates MUST always decrease costs.
 **/
class COINUTILSLIB_EXPORT CoinNodeHeap {
public:
  /**
   * Default Constructor.
   * Create the heap with space for nodes {0,...,numNodes-1}.
   * The cost of the nodes are set to infinity.
   **/
  explicit CoinNodeHeap(size_t numNodes);

  /**
   * Destructor
   **/
  ~CoinNodeHeap();

  /**
   * Reset the values in the heap.
   * All costs are set to infinity.
   **/
  void reset();

  /**
   * Update, always in decreasing order, the cost of a node.
   **/
  void update(size_t node, double cost);

  /**
   * Remove the next element, returning its cost.
   *
   * @param node used to store the element that was removed
   **/
  double removeFirst(size_t *node);

  /**
   * Check if the cost of the nodes are set to infinity.
   **/
  bool isEmpty() const;

private:
  /**
   * Priority queue itself
   **/
  std::vector<std::pair<size_t, double> > pq_;

  /**
   * Indicates the position of each node in pq
   **/
  std::vector<size_t> pos_;

  /**
   * Number of nodes of the heap.
   **/
  size_t numNodes_;
};


#endif //COINNODEHEAP_HPP
