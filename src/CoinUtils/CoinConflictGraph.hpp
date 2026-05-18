/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Abstract class for a Conflict Graph, see CoinStaticConflictGraph and
 * CoinDynamicConflictGraph for concrete implementations.
 *
 * @file CoinConflictGraph.hpp
 * @brief Abstract class for conflict graph
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef CONFLICTGRAPH_HPP
#define CONFLICTGRAPH_HPP

#include "CoinUtilsConfig.h"
#include <cstddef>
#include <vector>
#include <string>
#include <utility>

/**
 * Base class for Conflict Graph: a conflict graph
 * is a structure that stores conflicts between binary
 * variables. These conflicts can involve the original
 * problem variables or complementary variables.
 **/
class COINUTILSLIB_EXPORT CoinConflictGraph {
public:
  /**
   * Metadata describing the constraint responsible for an implied bound.
   */
  struct ImplicationRow {
    ImplicationRow() = default;
    ImplicationRow(std::string name, int idx) : rowName(std::move(name)), rowIndex(idx) {}
    std::string rowName;
    int rowIndex = -1;
  };

  /**
   * Captures a contradiction where different constraints imply opposite
   * fixes for the same binary variable.
   */
  struct BinaryBoundInfeasibility {
    size_t variableIndex = 0;
    std::string variableName;
    ImplicationRow fixedToZero;
    ImplicationRow fixedToOne;
  };

  CoinConflictGraph() { }

  /**
   * \brief Construct a conflict graph for a problem with the given number of columns.
   *
   * The conflict graph stores entries for each variable and its complement,
   * so the internal number of nodes will typically be `2 * _size`.
   * \param _size Number of columns in the mixed-integer linear program.
   */
  CoinConflictGraph(size_t _size);

  /**
   * \brief Copy-construct a conflict graph from another instance.
   * \param other Pointer to the conflict graph to be copied. If `nullptr`,
   * behavior is undefined.
   */
  CoinConflictGraph(const CoinConflictGraph *other);

  /**
   * Destructor
   **/
  virtual ~CoinConflictGraph();

  /**
   * \brief Check whether two nodes conflict.
   *
   * \param n1 First node index.
   * \param n2 Second node index.
   * \return `true` if there is an edge between \p n1 and \p n2 in the
   * conflict graph, `false` otherwise.
   */
  bool conflicting(size_t n1, size_t n2) const;

  /**
   * \brief Return all nodes conflicting with a given node.
   *
   * \param node Index of the node whose conflicts are requested.
   * \param temp Temporary storage area for storing conflicts; should have
   * space for all elements in the graph (see `size()`).
   * \param iv Auxiliary incidence array used to eliminate duplicates; it
   * should have length equal to the graph size and be initialized to `false`.
   * \return A pair `(n, ptr)` where `n` is the number of conflicting nodes
   * and `ptr` points to an array of conflicting node indices. The returned
   * array may be `temp` (if the temporary storage was used) or an internal
   * array owned by the conflict graph implementation.
   */
  std::pair< size_t, const size_t* > conflictingNodes ( size_t node, size_t* temp, char *iv ) const;

  /**
   * Density of the conflict graph:
   * (nConflicts / maxConflicts)
   **/
  double density() const;

  /**
   * Number of nodes in the conflict graph.
   **/
  size_t size() const;

  /**
   * Degree of a given node.
   **/
  virtual size_t degree( const size_t node ) const = 0;

  /**
   * Modified degree of a given node. The modified
   * degree of a node is the sum of its degree
   * with the degrees of its neighbors.
   **/
   virtual size_t modifiedDegree( const size_t node ) const = 0;

  /**
   * Minimum node degree.
   **/
  size_t minDegree( ) const;

  /**
   * Maximum node degree.
   **/
  size_t maxDegree( ) const;

  /**
   * Number of cliques stored explicitly.
   **/
  virtual size_t nCliques() const = 0;

  /**
   * Size of the i-th clique stored explicitly.
   **/
  virtual size_t cliqueSize( size_t idxClique ) const = 0;

  /**
   * Contents of the i-th clique stored explicitly.
   **/
  virtual const size_t *cliqueElements( size_t idxClique ) const = 0;

  /**
   * Return how many explicit cliques a node appears.
   **/
  virtual size_t nNodeCliques(size_t idxClique) const = 0;

  /**
   * Return which cliques a node appears.
   **/
  virtual const size_t *nodeCliques(size_t idxClique) const = 0;

  /**
   * Return the number of pairwise conflicts
   * stored for a node.
   **/
  virtual size_t nDirectConflicts( size_t idxNode ) const = 0;

  /**
   * List of pairwise conflicts (not stored as
   * cliques) for a node.
   **/
  virtual const size_t *directConflicts( size_t idxNode ) const = 0;

  /**
   * Recompute the degree of each node of the graph.
   **/
  void recomputeDegree();

  /**
   * Recompute the modified degree of each node
   * of the graph.
   **/
  void computeModifiedDegree();

  /**
   * Total number of conflicts stored directly.
   **/
  virtual size_t nTotalDirectConflicts() const = 0;

  /**
   * Total number of clique elements stored.
   **/
  virtual size_t nTotalCliqueElements() const = 0;

#ifdef CGRAPH_DEEP_DIVE
  /**
   * \brief Validate the conflict graph entries against solver data and a MIP start.
   *
   * This routine verifies that the entries stored in the conflict graph are
   * consistent with the provided solver column information and with a
   * candidate MIP start (partial assignment). It is intended for debugging
   * and deep-inspection builds only.
   *
   * \param numCols Number of columns in the model.
   * \param colTypes Array of length \p numCols indicating column types (expects 'B' for binary variables).
   * \param colLower Array of length \p numCols containing column lower bounds.
   * \param colUpper Array of length \p numCols containing column upper bounds.
   * \param colNames Vector of column names (size should be \p numCols).
   * \param mipStart Vector of (variable name, value) pairs representing a candidate MIP start expected to be conflict-free.
   */
  void validateConflictGraphUsingFeasibleSolution(
      size_t numCols,
      const char *colTypes,
      const double *colLower,
      const double *colUpper,
      const std::vector<std::string> &colNames,
      const std::vector< std::pair< std::string, double > > &mipStart) const;
#endif // CGRAPH_DEEP_DIVE

  /**
   * Print summarized information about
   * the conflict graph.
   **/
  void printSummary() const;

  /**
   * Return all infeasibilities detected while building the conflict graph.
   */
  const std::vector< BinaryBoundInfeasibility > &infeasibleImplications() const;

  /**
   * Set the the minimum size of a clique
   * to be explicitly stored as a clique
   * (not pairwise).
   **/
  static void setMinCliqueRow(size_t minClqRow);

  /**
   * Return the the minimum size of a clique
   * to be explicitly stored as a clique
   * (not pairwise).
   **/
  static size_t getMinCliqueRow();

  /** Set the maximum number of cliques to store during construction.
   *  Once reached, further large cliques are skipped (small cliques
   *  are still expanded into pairwise edges).  Default: 1200. */
  static void setMaxCliques(size_t maxClq);

  /** Return the maximum cliques limit. */
  static size_t getMaxCliques();

protected:
  /**
   * Parameter that controls the minimum size of
   * a clique to be explicitly stored as a clique
   * (not pairwise).
   **/
  static size_t minClqRow_;

  /** Maximum number of cliques to store. */
  static size_t maxCliques_;

  void registerBoundImplicationInfeasibility(const BinaryBoundInfeasibility &info);

  /**
   * Sets the degree of a node
   *
   * @param idxNode index of the node
   * @param deg degree of the node
   **/
  virtual void setDegree( size_t idxNode, size_t deg ) = 0;

  /**
   * Sets the modified degree of a node
   *
   * @param idxNode index of the node
   * @param mdegree modified degree of the node
   **/
  virtual void setModifiedDegree( size_t idxNode, size_t mdegree ) = 0;

  /**
   * Checks if two nodes are conflicting, considering
   * only the conflicts explicitly stored as cliques.
   **/
  bool conflictInCliques( size_t idxN1, size_t idxN2) const;

  /**
   * Initializes the structures of the conflict graph.
   **/
  void iniCoinConflictGraph(size_t _size);

  /**
   * Default constructor
   * @param other conflict graph to be copied
   **/
  void iniCoinConflictGraph(const CoinConflictGraph *other);

  /**
   * Number of nodes of the graph.
   **/
  size_t size_;

  /**
   * Number of conflicts (edges)
   * of the graph.
   **/
  size_t nConflicts_;

  /**
   * Maximum number of conflicts that
   * the graph can have.
   **/
  double maxConflicts_; //this number could be large, storing as double

  /**
   * Density of the graph
   **/
  double density_;

  /**
   * Indicates if the modified degree of the nodes
   * must be recomputed.
   **/
  bool updateMDegree;

  /**
   * Minimum degree of the nodes.
   **/
  size_t minDegree_;

  /**
   * Maximum degree of the nodes.
   **/
  size_t maxDegree_;

  std::vector< BinaryBoundInfeasibility > infeasibleImplications_;

#ifdef CGRAPH_STATS
public:
  struct RowTypeStats {
    size_t nRows = 0;
    double totalTime = 0.0;
    size_t rowsWithConflicts = 0;
    size_t rowsWithFixings = 0;
  };

  /** Per-row profile: stats grouped by (nz_bucket, sense, abs_rhs_bucket). */
  struct RowProfileStats {
    int nzBucket = 0;     // 0:<4, 1:<8, 2:<16, 3:<32, 4:<64, 5:<128, 6:<256, 7:>=256
    char sense = '?';     // 'L','G','E','R'
    int rhsBucket = 0;    // 0:0, 1:1, 2:2-5, 3:6-20, 4:21-100, 5:>100
    size_t nRows = 0;
    double totalTime = 0.0;
    size_t rowsWithConflicts = 0;
    size_t rowsWithFixings = 0;
    size_t totalConflictsFound = 0; // sum of clique sizes across rows
  };

  const RowTypeStats *rowTypeStats() const { return rowTypeStats_; }
  const std::vector<RowProfileStats> &rowProfileStats() const { return rowProfileStats_; }

protected:
  RowTypeStats rowTypeStats_[16]; // indexed by CoinRowType
  std::vector<RowProfileStats> rowProfileStats_;
#endif
};

#endif // CONFLICTGRAPH_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
