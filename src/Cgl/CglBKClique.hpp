/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for separating clique cuts.
 * It uses the Bron-Kerbosch (BK) algorithm.
 *
 * @file CglBKClique.hpp
 * @brief Clique cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef _CglBKClique_h_
#define _CglBKClique_h_

#include <CglCutGenerator.hpp>
#include "CoinBronKerbosch.hpp"
#include <atomic>

class CoinCliqueList;
class CoinConflictGraph;

/**
 * Class responsible for separating clique cuts.
 * It uses the Bron-Kerbosch algorithm.
 **/
class CGLLIB_EXPORT CglBKClique : public CglCutGenerator {
public:
  /**
   * Default constructor
   **/
  CglBKClique();

  /**
   * Copy constructor
   **/
  CglBKClique(const CglBKClique &rhs);

  /**
   * Clone
   **/
  virtual CglCutGenerator *clone() const;

  /**
   * Generate clique cuts for the model data contained
   * in si. The generated cuts are inserted into and returned
   * in the collection of cuts cs.
   **/
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info = CglTreeInfo());

  /**
   * Destructor
   **/
  virtual ~CglBKClique();

  /**
   * Refresh the conflict graph if necessary.
   **/
  virtual void refreshSolver(OsiSolverInterface *solver);

  /**
   * Set the maximum number of recursive calls to be made
   * by BK algorithm.
   **/
  void setMaxCallsBK(size_t maxCallsBK);

  /**
   * Return the maximum number of recursive calls to be made
   * by BK algorithm.
   **/
  size_t getMaxCallsBK() const { return maxCallsBK_; }

  /**
   * Set the wall-clock budget (in seconds) for a single generateCuts() call.
   * BK recursion is curtailed when the deadline approaches.
   * 0.0 = no limit (default).
   **/
  void setMaxSeconds(double maxSeconds) { maxSeconds_ = maxSeconds; }

  /**
   * Return the current wall-clock budget.
   **/
  double getMaxSeconds() const { return maxSeconds_; }

  /**
   * Set the maximum number of induced-subgraph vertices passed to BK.
   * When the LP solution produces more fractional/complement vertices than
   * this limit, they are ranked by LP-value * sqrt(conflict-degree) and only
   * the top maxInducedSize_ are kept.  This bounds the O(n^2) bitstring
   * construction cost in CoinBronKerbosch.
   * 0 = unlimited (default 10000).
   **/
  void setMaxInducedSize(size_t maxInducedSize) { maxInducedSize_ = maxInducedSize; }

  /**
   * Return the current induced-subgraph size cap.
   **/
  size_t getMaxInducedSize() const { return maxInducedSize_; }

  /**
   * Set the strategy used to extend cliques:
   * 0 = no extension; 1 = random; 2 = max degree;
   * 3 = max modified degree; 4 = reduced cost (inversely proportional);
   * 5 = reduced cost (inversely proportional) + modified degree.
   **/
  void setExtendingMethod(size_t extMethod);

  /**
   * Return the strategy used to extend cliques.
   **/
  size_t getExtendingMethod() const { return extMethod_; }

  /**
   * Return the number of calls made by BK algorithm.
   **/
  size_t getNumCallsBK() const { return callsBK_; }

  /**
   * Return the minimum value that a variable must have
   * in the LP relaxation to be considered in the cut
   * separation.
   **/
  const double getMinFrac() const { return minFrac_; }

  /**
   * Return the minimum violation of a clique to be
   * stored by the cut separator.
   **/
  const double getMinViol() const { return minViol_; }

  /**
   * Set the minimum value that a variable must have
   * in the LP relaxation to be considered in the cut
   * separation.
   **/
  void setMinFrac(const double minFrac);

  /**
   * Set the minimum violation of a clique to be
   * stored by the cut separator.
   **/
  void setMinViol(const double minViol);

  /**
   * Set the pivoting strategy used in BK algorithm
   **/
  void setPivotingStrategy(const CoinBronKerbosch::PivotingStrategy pivotingStrategy);

  /**
   * Number of cuts separated.
   **/
  static std::atomic< size_t > sepCuts_;

  /**
   * Execution time spent for the clique
   * cut separator.
   **/
  static std::atomic< double > sepTime_;

private:
  /**
   * Check if it is necessary realloc the memory
   * for the data structures.
   **/
  void checkMemory(const size_t newNumCols);

  /**
   * Execute the clique cut separation.
   * Return a list of violated cliques.
   **/
  CoinCliqueList *separateCliques(const OsiSolverInterface &si);

  /**
   * Execute the clique extension procedure.
   * Return a list of violated cliques.
   **/
  CoinCliqueList *extendCliques(const OsiSolverInterface &si, const CoinCliqueList *initialCliques);

  /**
   * Insert the violated cuts in OsiCuts cs.
   **/
  void insertCuts(const OsiSolverInterface &si, const CglTreeInfo &info, const CoinCliqueList *cliques, OsiCuts &cs);

  /**
   * Capacity of storage of the data structures.
   **/
  size_t cap_;

  /**
   * Reduced costs of the variables.
   **/
  double *rc_;

  /**
   * Auxiliary arrays used to store the indexes
   * of a clique.
   **/
  int *idxs_, *idxMap_;

  /**
   * Auxiliary array used to store the coefficients
   * of a clique.
   **/
  double *coefs_;

  /**
   * Auxiliary array that stores the indexes
   * to construct a induced subgraph.
   **/
  size_t *inducedVert_;

  /**
   * Auxiliary array for storing the current
   * clique.
   **/
  size_t *currClq_;

  /**
   * Weight of each vertex.
   **/
  double *vertexWeight_;

  /**
   * Minimum value that a variable must have
   * in the LP relaxation to be considered in the cut
   * separation.
   **/
  double minFrac_;

  /**
   * Minimum violation of a clique to be
   * stored by the cut separator.
   **/
  double minViol_;

  /**
   * Convert minViol_ to a integer value,
   * to made further computations easier.
   **/
  double minWeight_;

  /**
   * Pivoting strategy employed in the BK algorithm.
   **/
  CoinBronKerbosch::PivotingStrategy pivotingStrategy_;

  /**
   * Strategy used to extend cliques.
   **/
  size_t extMethod_;

  /**
   * Maximum number of recursive calls to be made
   * by BK algorithm.
   **/
  size_t maxCallsBK_;

  /**
   * Number of calls made by BK algorithm.
   **/
  size_t callsBK_;

  /**
   * Check if BK ran completely (without
   * stopping by maxCallsBK).
   **/
  bool completeBK_;

  /**
   * Wall-clock budget in seconds for a single generateCuts() call.
   * 0.0 = no limit.
   **/
  double maxSeconds_;

  /**
   * Cap on the number of induced-subgraph vertices sent to BK.
   * When the LP produces more fractional (original + complement) vertices,
   * they are ranked by LP-value * sqrt(conflict-degree) and only the top
   * maxInducedSize_ survive.  0 = unlimited.
   **/
  size_t maxInducedSize_;

  /**
   * Auxiliary structure used to temporary
   * store a cut.
   **/
  OsiRowCut osrc_;
};

#endif // CglBKClique_HPP
