/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for extending cliques
 *
 * @file CoinCliqueExtender.hpp
 * @brief Clique extender
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef COINCLIQUEEXTENDER_HPP
#define COINCLIQUEEXTENDER_HPP

#include <vector>
class CoinCliqueList;
class CoinConflictGraph;

#include "CoinUtilsConfig.h"
#include <cstddef>

/**
 * A class to extend cliques. It stores the cliques
 * extended in a clique list.
 **/
class COINUTILSLIB_EXPORT CoinCliqueExtender {
public:
  /**
   * Default constructor.
   *
   * @param cgraph conflict graph
   * @param extMethod method used to extend cliques (optional).
   * Values: 0 = no extension; 1 = random; 2 = max degree;
   * 3 = max modified degree; 4 = reduced cost (inversely proportional);
   * 5 = reduced cost (inversely proportional) + modified degree.
   * @param rc reduced cost (optional)
   **/
  explicit CoinCliqueExtender(const CoinConflictGraph *cgraph, size_t extMethod = 4, const double *rc = NULL,
          const double maxRC = 100.0);

  /**
   * Destructor
   **/
  ~CoinCliqueExtender();

  /**
   * Try to extend a clique. In case of success,
   * return true and store the clique extended
   * in a clique list.
   *
   * @param clqIdxs indexes of the clique to be extended
   * @param clqSize size of the clique to be extended
   **/
  bool extendClique(const size_t *clqIdxs, const size_t clqSize);

  /**
   * Return the number of cliques stored in the
   * list of cliques extended.
   **/
  size_t nCliques() const;

  /**
   * Return the i-th clique extended
   **/
  const size_t* getClique(const size_t i) const;

  /**
   * Return the size of the i-th clique extended
   **/
  size_t getCliqueSize(const size_t i) const;

  /**
   * Set the maximum number of candidates to enter
   * a clique.
   **/
  void setMaxCandidates(const size_t maxCandidates);

private:
  /**
   * Construct a list of candidates to enter the clique.
   *
   * @param clqIdxs indexes of the clique to be extended
   * @param clqSize size of the clique to be extended
   **/
  void fillCandidates(const size_t *clqIdxs, const size_t clqSize);

  /**
   * Randomly select the vertices to be inserted in the clique.
   *
   * @param clqIdxs indexes of the clique to be extended
   * @param clqIdxs size of the clique to be extended
   **/
  bool randomExtension(const size_t *clqIdxs, const size_t clqSize);

  /**
   * Select the vertices to be inserted in the clique using a greedy
   * strategy. Consider the method to extend cliques passed in
   * the constructor.
   *
   * @param clqIdxs indexes of the clique to be extended
   * @param clqIdxs size of the clique to be extended
   **/
  bool greedySelection(const size_t *clqIdxs, const size_t clqSize, const double *costs);

  /**
   * Conflict graph
   **/
  const CoinConflictGraph *cgraph_;

  /**
   * List of extended cliques
   **/
  CoinCliqueList *extendedCliques_; //stores the extended cliques

  /**Clique extension method: 0 = no extension;1 = random;
   * 2 = max degree; 3 = max modified degree;
   * 4 = reduced cost (inversely proportional);
   * 5 = reduced cost (inversely proportional) + modified degree.
  **/
  size_t extMethod_;

  /**
   * Maximum size of the candidates list
   **/
  size_t maxCandidates_;

  /**
   * Candidates list, used to temporarily store candidate vertices.
   **/
  std::vector<size_t> candidates_;
  size_t nCandidates_;

  /**
   * Auxiliary array used to temporarily
   * store a new extended clique.
   **/
  std::vector<size_t> newClique_;
  size_t nNewClique_;

  /**
   * Array containing the cost of each vertex,
   * computed according to the clique extension method.
   * Used to greedy select the candidates to enter a clique.
   **/
  std::vector<double> costs_;

  /**
   * Auxiliary incidence vectors
   **/
  std::vector<char> iv_, iv2_;

  /**
   * Array containing the reduced cost associated
   * with each vertex.
   **/
  const double *rc_;

  /**
   * Maximum value of the reduced cost that
   * a vertex must have to enter the candidates
   * list.
   **/
  double maxRC_;
};


#endif //COINCLIQUEEXTENDER_HPP
