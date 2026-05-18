/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for storing a pool of cuts, removing
 * the repeated and dominated ones. It also filters the cuts
 * according to their scores. Given a solution x* of the LP
 * relaxation of a MILP, the score of a cut C is calculated as
 * S(C) = viol(C) / actv(C), where viol(C) is the violation
 * of the cut with respect to x* and actv(C) is the
 * number of variables in C whose values x* are greater than
 * zero. An auxiliary array is used to identify, for each variable,
 * the cut in the pool with the best score that contains this variable.
 * A cut is only inserted into the cut pool if it has the best score for
 * at least one variable.
 *
 * @file CoinCutPool.hpp
 * @brief Class for storing a pool of cuts
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef COINCUTPOOL_HPP
#define COINCUTPOOL_HPP

#include "CoinUtilsConfig.h"
#include <cstddef>
#include <vector>

/**
 * Class for representing a cut.
 **/
class COINUTILSLIB_EXPORT CoinCut {
public:
  /**
   * Default constructor
   *
   * @param idxs indexes of variables of the cut
   * @param coefs coefficients of the variables the cut
   * @param nz size of the cut
   * @param rhs right-hand side of the cut
   **/
  CoinCut(const int *idxs, const double *coefs, size_t nz, double rhs);

  /**
   * Destructor
   **/
  ~CoinCut();

  /**
   * Return the indexes of the variables of the cut.
   **/
  const int* idxs() const { return idxs_.data(); }

  /**
   * Return the coefficients of the variables of the cut.
   **/
  const double* coefs() const { return  coefs_.data(); }

  /**
   * Return the size of the cut.
   **/
  int size() const { return idxs_.size(); }

  /**
   * Return the right-hand side of the cut.
   **/
  double rhs() const { return rhs_; }

  /**
   * Check if the cut dominates another one.
   *
   * @param other cut to be checked.
   **/
  bool dominates(const CoinCut *other) const;

private:
  /**
   * indexes of the variables of the cut
   **/
  std::vector<int> idxs_;

  /**
   * coefficients of the variables of the cut
   **/
  std::vector<double> coefs_;

  /**
   * right-hand side of the cut
   **/
  double rhs_;
};

/**
 * Class for storing a pool of cuts, removing
 * the repeated and dominated ones. It also filters the cuts
 * according to their scores.
 **/
class COINUTILSLIB_EXPORT CoinCutPool {
public:
  /**
   * Default constructor
   *
   * @param x current solution of the LP relaxation of the MILP.
   * @param numCols number of variables of the MILP.
   **/
  CoinCutPool(const double *x, int numCols);

  /**
   * Destructor
   **/
  ~CoinCutPool();

  /**
   * Return the number of cuts in the pool.
   **/
  size_t numCuts() const;

  /**
   * Return the indexes of the variables of
   * the i-th cut in the pool.
   **/
  const int* cutIdxs(size_t i) const;

  /**
   * Return the coefficients of the variables
   * of the i-th cut in the pool.
   **/
  const double* cutCoefs(size_t i) const;

  /**
   * Return the size of the i-th cut in the
   * pool.
   **/
  int cutSize(size_t i) const;

  /**
   * Return the right-hand side of the
   * i-th cut in the pool.
   **/
  double cutRHS(size_t i) const;

  /**
   * Try to add a cut in the pool. Return true
   * if the cut was added.
   *
   * @param idxs indexes of the variables of the
   * cut to be added.
   * @param coefs coefficients of the variables of
   * the cut to be added.
   * @param nz size of the cut to be added
   * @param rhs right-hand side of the cut to be added
   **/
  bool add(const int *idxs, const double *coefs, int nz, double rhs);

  /**
   * Remove dominated cuts.
   **/
  void removeDominated();

  /**
   * Remove cuts that are previous deleted.
   * Must be called before use the cuts in
   * the pool.
   **/
  void removeNullCuts();

private:
  /**
   * Compute the number of variables for which the
   * cut has the best score. Update the structures
   * to reflect this calculation. Return the number
   * of variables for which the cut has the best score.
   **/
    size_t updateCutFrequency(const CoinCut *cut);

  /**
   * Compute the score of the cut.
   **/
  double calculateFitness(const CoinCut *cut) const;

  /**
   * Check if it is necessary expand the memory
   * allocated for storing the cuts.
   **/
  void checkMemory();

  /**
   * Check the dominance relation between two cuts
   * in the pool. Return 0 if cut idxA dominates
   * cut idxB, 1 if cut idxB dominates cut idxA,
   * or 2 if none is not dominated by the other.
   **/
  int checkCutDomination(size_t idxA, size_t idxB);

  /**
    * Array of pointers to the cuts stored in the pool.
    **/
  std::vector<CoinCut *> cuts_;

  /**
   * Number of cuts stored in the pool.
   **/
  size_t nCuts_;

  /**
   * Current capacity to store cuts in the pool.
   **/
  size_t cutsCap_;

  /**
   * For each cut, stores the number of variables in
   * which it has the best score.
   **/
  std::vector<size_t> cutFrequency_;

  /**
   * Score of each cut
   **/
  std::vector<double> cutFitness_;

  /**
   * For each variable, stores the index of the cut
   * with the best score that contains this variable.
   **/
  std::vector<int> bestCutByCol_;

  /**
   * Number of cuts that were deleted from the pool.
   **/
  size_t nullCuts_;

  /**
   * Current solution of the LP relaxation of the MILP.
   **/
  const double *x_;
};


#endif //COINCUTPOOL_HPP
