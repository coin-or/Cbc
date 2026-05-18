/**
 *
 * This file is part of the COIN-OR CoinUtils project
 *
 * CoinConflictGraph implementation which supports modifications.
 *
 * For a static conflict graph implemenation with faster queries
 * check CoinStaticConflictGraph.
 *
 * @file CoinDynamicConflictGraph.hpp
 * @brief CoinConflictGraph implementation which supports modifications.
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef DYNAMICCONFLICTGRAPH_H
#define DYNAMICCONFLICTGRAPH_H

#include <vector>
#include <map>
#include <utility>
#include <string>

#include "CoinUtilsConfig.h"
#include "CoinTypes.h"
#include "CoinConflictGraph.hpp"
#include "CoinAdjacencyVector.hpp"
#include "CoinTerm.hpp"
#include "CoinKnapsackRow.hpp"

class CoinPackedMatrix;
class CoinAdjacencyVector;
class CoinCliqueList;

/**
 * This a a conflict graph where conflicts can be added on the fly,
 * not optimized for memory usage.
 **/
class COINUTILSLIB_EXPORT CoinDynamicConflictGraph : public CoinConflictGraph
{
public:
  /**
   * Default constructor
   *
   * @param _size number of vertices of the
   * conflict graph
   **/
  CoinDynamicConflictGraph ( size_t _size );

  /**
   * Builds a conflict graph directly from a MILP instance.
   *
   * The constructor scans every constraint, complements variables as
   * needed, and detects both explicit cliques (stored in `largeClqs` or as
   * direct conflicts) and partial cliques that must be revisited later
   * through `cliqueDetection`. During this pass it also:
   *   - discounts fixed / continuous columns from each row so that only
   *     binary columns remain candidates for conflicts,
   *   - inserts trivial variable–complement conflicts,
   *   - stores temporary rows for further processing of small cliques, and
   *   - determines tighter bound recommendations that can later be queried
   *     via `updatedBounds()`.
   *
   * @param numCols number of original columns
   * @param colType column types (used to filter binary, semi-continuous,
   *                etc.)
   * @param colLB column lower bounds
   * @param colUB column upper bounds
   * @param matrixByRow row-wise constraint matrix
   * @param sense row senses for each row (G/E/L/R/N)
   * @param rowRHS right-hand side for each row
   * @param rowRange row ranges (only used for range rows)
   * @param primalTolerance numerical tolerance used while detecting nearly
   *                        tight inequalities (e.g., when testing a_j > b)
   * @param infinity threshold to interpret bounds as practically unbounded
   *                 when discounting continuous variables
   * @param colNames column names, used exclusively for diagnostic messages
   *                 while reporting newly inferred bounds
   **/
  CoinDynamicConflictGraph(
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
    const std::vector< std::string > &rowNames
  );

  /**
   * Destructor
   **/
  virtual ~CoinDynamicConflictGraph();

  /**
   * Add conflicts to a node to be stored directly (not as cliques).
   *
   * @param node index
   * @param nodeConflicts conflicts to be added
   * @param nConflicts number of conflicts to be added
   **/
  void addNodeConflicts( const size_t node, const size_t nodeConflicts[], const size_t nConflicts );

  /**
   * Add a clique (it will be stored explicitly as a clique).
   **/
  void addClique( size_t size, const size_t elements[] );

  /**
   * Return the number of cliques stored explicitly.
   **/
  size_t nCliques() const;

  /**
   * Return the contents of the i-th clique stored explicitly.
   **/
  const size_t* cliqueElements(size_t idxClique) const ;

  /**
   * Return the size of the i-th clique stored explicitly.
   **/
  size_t cliqueSize( size_t idxClique ) const;

  /**
   * Return how many explicit cliques a node appears.
   **/
  size_t nNodeCliques(size_t idxNode) const;

  /**
   * Return which cliques a node appears.
   **/
  const size_t *nodeCliques(size_t idxNode) const;

  /**
   * Return the degree of a given node.
   **/
  size_t degree( const size_t node ) const;

  /**
   * Return the modified degree of a given node.
   **/
  size_t modifiedDegree( const size_t node ) const;

  /**
   * Total number of clique elements stored.
   **/
  size_t nTotalCliqueElements() const;

  /**
   * Return the number of pairwise conflicts
   * stored for a node.
   **/
  size_t nDirectConflicts( size_t idxNode ) const;

  /**
   * List of pairwise conflicts (not stored as
   * cliques) for a node.
   **/
  const size_t *directConflicts( size_t idxNode ) const;

  /**
   * Move the direct-conflict vector for a node (transfers ownership).
   **/
  std::vector<size_t> moveDirectConflicts(size_t idxNode);

  /**
   * Move a clique vector (transfers ownership).
   **/
  std::vector<size_t> moveClique(size_t idxClique);

  /**
   * Total number of conflicts stored directly.
   **/
  virtual size_t nTotalDirectConflicts() const;

  /**
   * Add a clique as pairwise conflicts.
   **/
  void addCliqueAsNormalConflicts( const size_t idxs[], const size_t len );

  /**
   * Insert trivial conflicts between each binary variable and its
   * complemented counterpart.
   */
  void addVariableComplementConflicts(int numCols, const char *colType);

  /**
   * Recommended tighter bounds for some variables
   *
   * The construction of the conflict graph may discover new tighter
   * bounds for some variables.
   *
   * @return a vector of updated bounds with the format (idx, (lb, ub))
   **/
  const std::vector< std::pair< size_t, std::pair< double, double > > > &updatedBounds();

  /**
   * Print information about the conflict graph.
   **/
  void printInfo() const;

private:
  /**
   * Sets the degree of a node.
   *
   * @param idxNode index of the node
   * @param deg degree of the node
   **/
  void setDegree( size_t idxNode, size_t deg );

  /**
   * Sets the modified degree of a node.
   *
   * @param idxNode index of the node
   * @param deg degree of the node
   **/
  void setModifiedDegree(size_t idxNode, size_t mdegree);

  /**
   * Return the size of the i-th clique
   **/
  size_t getClqSize( size_t idxClique ) const;

  /**
   * Detect all clique conflicts implied by a single processed constraint.
   *
   * The vector `columns` must contain only binary (possibly complemented)
   * columns that survived the preprocessing pass and be sorted in
   * nondecreasing order of their nonnegative coefficients. Given the
   * adjusted right-hand side `rhs`, the routine identifies the tightest
   * suffix whose pairwise sums violate the row, registers it as a clique via
   * `processClique`, and then iteratively seeds new cliques by combining
   * earlier columns with the suffix found through a bounded binary search.
   *
   * @param columns row entries (column index, coefficient) already filtered
   *                and sorted for clique detection
   * @param nz number of relevant columns stored in `columns`
   * @param rhs adjusted right-hand side after removing continuous/fixed parts;
   *            callers typically add the primal tolerance before invoking this
   *            routine so that comparisons treat nearly tight inequalities as
   *            conflicting.
   **/
  void cliqueDetection(const std::vector<CoinTerm> &columns, size_t nz, const double rhs);

  /**
   * Add a clique. It will be stored explicitly or not
   * depending on its size.
   **/
  void processClique( const size_t idxs[], const size_t size );

  /**
   * Promote a node's "small" cliques (kept temporarily in `smallCliques`)
   * to explicit pairwise conflicts.
   *
   * The routine receives the index of the node, the list (`scn`) with the
   * identifiers of the small cliques that contain that node, and uses the
   * scratch array `iv` to flag neighbors already implied either by previous
   * direct conflicts or by large cliques.  It then scans every referenced
   * small clique and inserts only the missing neighbors via
   * `fastAddNeighbor`, ensuring duplicates are not created and finally
   * restoring the scratch markers.
   *
   * @param node node whose small cliques are being materialized
   * @param scn array with the indices of the node's small cliques inside
   *           `smallCliques`
   * @param nscn number of entries in `scn`
   * @param smallCliques container that temporarily stores all detected small
   *                     cliques
   * @param iv scratch boolean array (size equal to the graph order) used to
   *           mark which neighbors are already present while processing
   **/
  void processSmallCliquesNode(
    size_t node,
    const size_t scn[],
    const size_t nscn,
    const CoinCliqueList *smallCliques,
    char *iv );

  /**
   * Cache a constraint for deferred clique detection.
   *
   * During the preprocessing pass the constructor separates rows that are not
   * immediate cliques (typically because only a suffix violates the rhs) and
   * stores the surviving binary columns together with the adjusted rhs in a
   * temporary buffer.  Later, `cliqueDetection()` iterates over this buffer to
   * mine all remaining implied cliques.  The method copies the first `nz`
   * `(column, coefficient)` pairs from `columns` to `tRowElements` and records
   * the corresponding right-hand side in `tRowRHS`.
   *
   * @param nz number of filtered binary entries that should be cached
   * @param columns workspace holding the filtered `(column, coefficient)`
   *                pairs; only the first `nz` entries are consumed and copied
   * @param rhs adjusted right-hand side after removing fixed/continuous
   *            contributions; stored verbatim for later reuse
   */
  void addTmpRow(
    size_t nz,
    const CoinTerm els[],
    double rhs);


  /**
   * Registers that a column must be fixed to either 0 or 1 in
   * updated bounds list (depending on whether it is an original
   * or complemented column).
   */
  void addFixingBound(size_t columnIdx,
                      size_t numCols,
                      size_t idxRow,
                      const std::vector<std::string> &colNames,
                      const std::vector<std::string> &rowNames
                    );

  /**
   * Extracts all fixed variables detected while processing a knapsack row
   * and registers their implied bounds.
   */
  void processKnapsackRowFixings(
    const CoinKnapsackRow &knapsackRowHelper,
    size_t numCols,
    size_t idxRow,
    const std::vector<std::string> &colNames,
    const std::vector<std::string> &rowNames);

#ifdef CGRAPH_DEEP_DIVE
  /**
   * Helper used for debugging a specific row.
   */
  void debugRowDetails(
    size_t idxRow,
    const char* sense,
    const double* rowRHS,
    const double* rowRange,
    const int *idxs,
    const double *coefs,
    const CoinBigIndex *start,
    const int *length,
    const std::vector<std::string> &colNames,
    const char* colType,
    const double* colLB,
    const double* colUB) const;
#endif // CGRAPH_DEEP_DIVE

  /**
   * Conflicts stored directly (not as cliques)
   **/
  CoinAdjacencyVector *conflicts;

  /**
   * Degree of the nodes
   **/
  std::vector<size_t> degree_;

  /**
   * Modified degree of the nodes
   **/
  std::vector<size_t> modifiedDegree_;

  /**
   * Conflicts stored as cliques
   **/
  CoinCliqueList *largeClqs;

  /**
   * Recommended tighter bounds for some variables.
   **/
  std::vector< std::pair< size_t, std::pair< double, double > > > newBounds_;

  /**
   * Temporary reservoir that holds the “small” cliques (those with fewer than
   * `CoinConflictGraph::minClqRow_` elements) detected during the initial row
   * scan.  They are buffered here so that, once all rows are processed and the
   * large cliques have been indexed, each node can revisit its small cliques
   * and materialize only the still-missing pairwise conflicts via
   * `processSmallCliquesNode()`.
   **/
  CoinCliqueList *smallCliques;

  /**
   * Rows that must be revisited later (stores the filtered column/value pairs).
   */
  std::vector<std::vector<CoinTerm> > tRowElements;

  /**
   * Tightened right-hand side associated with each temporary row above.
   */
  std::vector<double> tRowRHS;

  /**
   * Implications detected during conflict graph construction per column.
   *
   * For a column j, it stores implications (which can be j for original variable
   * or j + numCols for complemented variable) and which row caused the
   * implication (name and index).
   */
  std::map<int, std::map<int, std::pair<std::string, int>>> columnRowImplications_;
};

#endif // DYNAMICCONFLICTGRAPH_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
