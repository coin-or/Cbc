/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * CoinConflictGraph implementation which supports fast queries
 * but doesn't supports modifications.
 *
 * @file CoinStaticConflictGraph.hpp
 * @brief static CoinConflictGraph implementation with fast queries
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef STATICCONFLICTGRAPH_HPP
#define STATICCONFLICTGRAPH_HPP

#include "CoinUtilsConfig.h"
#include "CoinConflictGraph.hpp"
#include "CoinDynamicConflictGraph.hpp"

#include <algorithm>
#include <vector>

/**
 * Static conflict graph, optimized for memory usage and query speed,
 * not modifiable.
 **/
class COINUTILSLIB_EXPORT CoinStaticConflictGraph : public CoinConflictGraph
{
public:
  /**
   * Default constructor
   **/
  CoinStaticConflictGraph ( const CoinConflictGraph *cgraph );

  /**
   * Builds a static conflict graph directly from a MILP instance.
   *
   * Internally this constructor instantiates a `CoinDynamicConflictGraph`,
   * lets it scan every constraint to detect cliques, pairwise conflicts, and
   * tighter bounds, and then copies the resulting structure into the compact
   * static layout optimized for query speed. The tighter bounds discovered
   * during this pass can later be accessed through `updatedBounds()`.
   *
   * @param numCols number of original columns.
   * @param colType column types (used to filter the binary/semicontinuous
   *                variables that may participate in conflicts).
   * @param colLB column lower bounds.
   * @param colUB column upper bounds.
   * @param matrixByRow row-wise constraint matrix.
   * @param sense row senses for each constraint (G/E/L/R/N).
   * @param rowRHS right-hand sides for each row.
   * @param rowRange row ranges (only relevant for range rows).
   * @param primalTolerance feasibility tolerance applied while testing a
   *                        row for conflict generation.
   * @param infinity threshold used to treat bounds as practically unbounded
   *                 when discounting continuous columns.
   * @param colNames column names, only used for diagnostic messages when
   *                 reporting newly inferred bounds.
   **/
  CoinStaticConflictGraph(
          const int numCols,
          const char* colType,
          const double* colLB,
          const double* colUB,
          const CoinPackedMatrix* matrixByRow,
          const char* sense,
          const double* rowRHS,
          const double* rowRange,
          const double primalTolerance,
          const double infinity,
          const std::vector<std::string> &colNames,
          const std::vector<std::string> &rowNames
        );

  /**
   * Clone a conflict graph.
   **/
  CoinStaticConflictGraph *clone() const;

  /**
   * Constructor to create an induced subgraph
   *
   * @param cgraph conflict graph
   * @param n number of elements in the induced subgraph
   * @param elements indexes of nodes in the induced subgraph
   **/
  CoinStaticConflictGraph( const CoinConflictGraph *cgraph, const size_t n, const size_t elements[] );

  /**
   * Return the number of pairwise conflicts
   * stored for a node.
   **/
  virtual size_t nDirectConflicts( size_t idxNode ) const;

  /**
   * List of pairwise conflicts (not stored as
   * cliques) for a node.
   **/
  virtual const size_t *directConflicts( size_t idxNode ) const;

  /**
   * Return the number of cliques stored explicitly.
   **/
  virtual size_t nCliques() const;

  /**
   * Return the contents of the i-th clique stored explicitly.
   **/
  virtual const size_t *cliqueElements( size_t idxClique ) const;

  /**
   * Return the size of the i-th clique stored explicitly.
   **/
  virtual size_t cliqueSize( size_t idxClique ) const;

  /**
   * Return how many explicit cliques a node appears.
   **/
  size_t nNodeCliques( size_t idxNode ) const;

  /**
   * Return which cliques a node appears.
   **/
  const size_t *nodeCliques( size_t idxNode ) const;

  /**
   * degree of a given node
   */
  virtual size_t degree( const size_t node ) const;

  /**
   * Return the modified degree of a given node.
   **/
  virtual size_t modifiedDegree( const size_t node ) const;

  /**
   * Total number of conflicts stored directly.
   **/
  virtual size_t nTotalDirectConflicts() const;

  /**
   * Total number of clique elements stored.
   **/
  virtual size_t nTotalCliqueElements() const;

  /**
   * Destructor
   **/
  virtual ~CoinStaticConflictGraph();

  /**
   * Fast non-virtual conflict check.
   * Hides the base-class method for callers that use
   * CoinStaticConflictGraph directly (avoids virtual dispatch).
   **/
  bool conflicting(size_t n1, size_t n2) const {
    if (n1 == n2)
      return false;

    // Check direct conflicts (sorted arrays → binary search)
    const auto &dc1 = conflicts_[n1];
    const auto &dc2 = conflicts_[n2];
    if (dc1.size() < dc2.size()) {
      if (std::binary_search(dc1.begin(), dc1.end(), n2))
        return true;
    } else {
      if (std::binary_search(dc2.begin(), dc2.end(), n1))
        return true;
    }

    // Check clique conflicts
    const auto &nc1 = nodeCliques_[n1];
    const auto &nc2 = nodeCliques_[n2];
    size_t nodeToSearch;
    const std::vector<size_t> *clqList;
    if (nc1.size() < nc2.size()) {
      clqList = &nc1;
      nodeToSearch = n2;
    } else {
      clqList = &nc2;
      nodeToSearch = n1;
    }
    for (size_t idxClq : *clqList) {
      const auto &clq = cliques_[idxClq];
      if (std::binary_search(clq.begin(), clq.end(), nodeToSearch))
        return true;
    }
    return false;
  }

  /**
   * Recommended tighter bounds for some variables
   *
   * The construction of the conflict graph may discover new tighter
   * bounds for some variables.
   *
   * @return a vector of updated bounds with the format (idx, (lb, ub))
   **/
  const std::vector< std::pair< size_t, std::pair< double, double > > > &updatedBounds() const;

private:
  /**
   * Recommended tighter bounds for some variables.
   **/
  std::vector< std::pair< size_t, std::pair< double, double > > > newBounds_;

  /**
   * Check if a clique contains a node.
   **/
  bool nodeInClique( size_t idxClique, size_t node ) const;

  /**
   * Sets the degree of a node.
   *
   * @param idxNode index of the node
   * @param deg degree of the node
   **/
  virtual void setDegree(size_t idxNode, size_t deg);

  /**
   * Sets the modified degree of a node.
   *
   * @param idxNode index of the node
   * @param deg degree of the node
   **/
  virtual void setModifiedDegree(size_t idxNode, size_t mdegree);

  /**
   * Initializes the structures of the conflict graph.
   **/
  void iniCoinStaticConflictGraph(const CoinConflictGraph *cgraph);

  /**
   * Number of pairwise conflicts stored.
   **/
  size_t nDirectConflicts_;

  /**
   * Number of elements considering all cliques
   * stored explicitly.
   **/
  size_t totalCliqueElements_;

  /**
   * Degree of each node.
   **/
  std::vector<size_t> degree_;

  /**
   * Modified degree of each node.
   **/
  std::vector<size_t> modifiedDegree_;

  /**
   * Direct conflicts
   **/
  std::vector<std::vector<size_t> > conflicts_;

  /**
   * Array containing, for each node,
   * the cliques that contain this node.
   **/
  std::vector<std::vector<size_t> > nodeCliques_;

  /**
   * Elements of the cliques stored explicitly.
   **/
  std::vector<std::vector<size_t> > cliques_;
};

#endif // STATICCONFLICTGRAPH_H
