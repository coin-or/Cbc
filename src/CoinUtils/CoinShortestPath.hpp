/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class that implements Dijkstra's algorithm for
 * finding the shortest paths between nodes in a graph.
 * Specialized for sparse graphs.
 *
 * @file CoinShortestPath.hpp
 * @brief Shortest path algorithm
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef COINSHORTESTPATH_HPP
#define COINSHORTESTPATH_HPP

#include "CoinUtilsConfig.h"
#include <cstddef>
#include <utility>
#include <vector>

class CoinNodeHeap;

/**
 * Class that implements Dijkstra's algorithm for
 * finding the shortest paths between nodes in a graph.
 **/
class COINUTILSLIB_EXPORT CoinShortestPath {
public:
  /**
   * Default constructor
   *
   * @param nodes number of nodes in the graph
   * @param arcs number of arcs in the graph
   * @param arcStart array containing the start
   * position of the neighbors of each node
   * @param toNode array containing the neighbors
   * of all nodes
   * @param dist array containing the distance
   * between each node and its neighbors.
   **/
  CoinShortestPath(size_t nodes, size_t arcs, const size_t *arcStart, const size_t *toNode, const double *dist);

  /**
   * Destructor
   **/
  ~CoinShortestPath();

  /**
   * Execute the shortest path finder using Dijkstra's algorithm.
   * Find the shortest path from node origin to all the other
   * nodes.
   **/
  void find(const size_t origin);

  /**
   * Execute the shortest path finder using the Dijkstra algorithm.
   * Find the shortest path from node origin to node destination.
   **/
  void find(const size_t origin, const size_t destination);

  /**
   * Fill array indexes with all previous nodes which should be
   * steped to arrive at a given node (this node is not included).
   * Return how many nodes were filled in indexes.
   **/
  size_t path(size_t toNode, size_t *indexes);

  /**
   * Return the number of nodes in the graph.
   **/
  size_t numNodes() const { return nodes_; }

  /**
   * Return the number of arcs in the graph.
   **/
  size_t numArcs() const { return neighs_.size(); }

  /**
   * Return the length of the shortest path
   * from the origin node (passed in method find)
   * to node.
   **/
  double distance(size_t node) const;

  /**
   * Return the previous node of a given node
   * in the shortest path.
   **/
  size_t previous(size_t node) const;

  /**
   * Return an array containing the previous
   * nodes which should be steped to arrive
   * at a given node.
   **/
  const size_t* previous() const { return previous_.data(); }

private:
  /**
   * Number of nodes in the graph.
   **/
  size_t nodes_;

  /**
   * Array containing all neighbors
   **/
  std::vector<std::vector<std::pair<size_t, double> > > neighs_;

  /**
   * Length of the shortest paths.
   **/
  std::vector<double> dist_;

  /**
   * Array indexes with all previous nodes which
   * should be steped to arrive at a given node.
   **/
  std::vector<size_t> previous_;

  /**
   * Temporary storage for the shortest paths.
   **/
  std::vector<size_t> path_;

  /**
   * Monotone heap used in Dijkstra's algorithm.
   **/
  CoinNodeHeap *nh_;
};


#endif //COINSHORTESTPATH_HPP
