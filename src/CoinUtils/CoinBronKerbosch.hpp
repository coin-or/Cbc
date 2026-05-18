/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for finding cliques with weights
 * greater than a certain threshold. It implements the
 * Bron-Kerbosch Algorithm.
 *
 * @file CoinBronKerbosch.hpp
 * @brief Bron-Kerbosch Algorithm
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef COINBRONKERBOSCH_HPP
#define COINBRONKERBOSCH_HPP

#include <cstddef>
#include <vector>
#include "CoinUtilsConfig.h"

class CoinConflictGraph;
class CoinCliqueList;

/**
 * Struct used to represent a weighted-vertex.
 **/
struct BKVertex {
  /**
   * Index of the vertex
   **/
  size_t idx;

  /**
   * Weight of the vertex
   **/
  double weight;

  /**
   * Degree of the vertex
   **/
  size_t degree;

  /**
   * Fitness of the vertex. Used to select
   * the pivot vertex in BK.
   **/
  double fitness;
};

/**
 * Class responsible for finding cliques with weights
 * greater than a certain threshold. It implements the
 * Bron-Kerbosch Algorithm.
 **/
class COINUTILSLIB_EXPORT CoinBronKerbosch {
public:
  /**
   * Pivoting strategies used in BK algorithm.
   **/
  enum PivotingStrategy {
    Off = 0,
    Random = 1,
    Degree = 2,
    Weight = 3,
    ModifiedDegree = 4,
    ModifiedWeight = 5,
    ModifiedDegreeWeight = 6,
  };

  /**
   * Default constructor.
   *
   * @param cgraph conflict graph
   * @param weights array containing the weights for each vertex
   * @param pivotingStrategy pivoting strategy used in BK algorithm.
   **/
  CoinBronKerbosch(const CoinConflictGraph *cgraph, const double *weights,
    PivotingStrategy pivotingStrategy = PivotingStrategy::Weight);

  /**
   * Destructor
   **/
  ~CoinBronKerbosch();

  /**
   * Perform BK algorithm for finding cliques with weights
   * greater than a certain threshold.
   **/
  void findCliques();

  /**
   * Return the number of cliques found by BK algorithm.
   **/
  size_t nCliques() const;

  /**
   * Return the i-th clique found by BK algorithm.
   **/
  const size_t* getClique(const size_t i) const;

  /**
   * Return the size of the i-th clique found by BK algorithm.
   **/
  size_t getCliqueSize(const size_t i) const;

  /**
   * Return the weight of the i-th clique found by BK algorithm.
   **/
  double getCliqueWeight(const size_t i) const;

  /**
   * Set the minimum weight that a clique must have
   * to be detected and stored by BK algorithm.
   **/
  void setMinWeight(double minWeight);

  /**
   * Set the maximum number of recursive calls of BK algorithm.
   **/
  void setMaxCalls(size_t maxCalls);

  /**
   * Return true if BK algorithm ran completely, without
   * stopping by the maximum number of recursive calls.
   **/
  bool completedSearch() const;

  /**
   * Return the number of calls made by BK algorithm
   **/
  size_t numCalls() const;

private:
    /**
     * Compute the fitness of all vertices.
     *
     * @param array containing the weights for each vertex
     **/
    void computeFitness(const double *weights);

    /**
     * Return the sum of the weights of the vertices in P.
     * Store in u the vertice in P with the highest fitness.
     **/
    double weightP(size_t depth, size_t &u);

    /**
     * Call the BK algorithm
     **/
    void bronKerbosch(size_t depth);

    /**
     * Conflict graph
     **/
    const CoinConflictGraph *cgraph_;

    /**
     * Array of vertices
     **/
    std::vector<BKVertex> vertices_;

    /**
     * Size of each bit vector
     **/
    size_t sizeBitVector_;

    /**
     * Bit mask
     **/
    std::vector<size_t> mask_;

    /**
     * Conflict graph and its complement represented by bit vectors
     **/
    std::vector<std::vector<size_t> > cgBitstring_, ccgBitstring_;

    /**
     * Bitstring with all vertices
     **/
    std::vector<size_t> allIn_;

    /**
     * Set C of vertices that are part of the current clique
     **/
    std::vector<size_t> C_;
    size_t nC_;

    /**
     * Sum of the weights of the vertices in C
     **/
    double weightC_;

    /**
     * Set of the candidate vertices to enter in C
     * Contain a vector for each depth of BK algorithm.
     **/
    std::vector<std::vector<size_t> > P_;
    std::vector<size_t> nP_;

    /**
     * Set of all vertices that have already been considered
     * in earlier steps by BK. Contain a vector for each
     * depth of BK algorithm.
     **/
    std::vector<std::vector<size_t> > S_;
    std::vector<size_t> nS_;

    /**
     * Set of vertices P excluding the neighbors of a vertex v.
     * Contain a vector for each depth of BK algorithm.
     **/
    std::vector<std::vector<size_t> > L_;

    /**
     * Array that stores the weights of each clique found by BK
     **/
    std::vector<double> clqWeight_;

    /**
     * Current capacity of array clqWeight_
     **/
    size_t clqWeightCap_;

    /**
     * Minimum weight that a clique must have to be
     * detected and stored by BK algorithm.
     **/
    double minWeight_;

    /**
     * Number of calls made by BK algorithm
     **/
    size_t calls_;

    /**
     * Maximum number of recursive calls of BK algorithm.
     **/
    size_t maxCalls_;

    /**
     * Pivoting strategy used in BK algorithm.
     **/
    PivotingStrategy pivotingStrategy_;

    /**
     * If BK algorithm ran completely, without stopping
     * by the maximum number of recursive calls.
     **/
    bool completeSearch_;

    /**
     * Cliques found by BK algorithm
     **/
    CoinCliqueList *cliques_;
};


#endif //COINBRONKERBOSCH_HPP
