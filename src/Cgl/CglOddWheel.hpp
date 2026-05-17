/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 *
 * @file CglOddWheel.hpp
 * @brief Odd-wheel cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef _CglOddWheel_h_
#define _CglOddWheel_h_

#include "CglCutGenerator.hpp"
#include "CoinConflictGraph.hpp"

class CoinConflictGraph;

/**
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 **/
class CGLLIB_EXPORT CglOddWheel : public CglCutGenerator
{
public:
  /**
   * Number of cuts separated.
   **/
  static size_t sepCuts;

  /**
   * Execution time spent for the clique
   * cut separator.
   **/
  static double sepTime;

  /**
   * Default constructor
   *
   * @param extMethod strategy that will be used to lift odd cycles,
   * transforming them into odd wheels: 0 = no lifting, 1 = only one
   * variable as wheel center, 2 = a clique as wheel center.
   **/
  CglOddWheel(size_t extMethod = 2);

  /**
   * Copy constructor
   **/
  CglOddWheel(const CglOddWheel& rhs);

  /**
   * Clone
   **/
  virtual CglCutGenerator * clone() const;

  /**
   * Generate clique cuts for the model data contained
   * in si. The generated cuts are inserted into and returned
   * in the collection of cuts cs.
   **/
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info = CglTreeInfo() );

  /**
   * Destructor
   **/
  virtual ~CglOddWheel();

  /**
   * Refresh the conflict graph if necessary.
   **/
  virtual void refreshSolver(OsiSolverInterface *solver);

  /**
   * Set the strategy that will be used to lift odd cycles,
   * transforming them into odd wheels: 0 = no lifting, 1 = only one
   * variable as wheel center, 2 = a clique as wheel center.
   **/
  void setExtendingMethod(size_t extMethod);

  /**
   * Return the strategy used to lift odd cycles.
   **/
  size_t getExtendingMethod() const { return extMethod_; }

private:
  /**
   * Check if it is necessary realloc the memory
   * for the data structures.
   **/
  void checkMemory(const size_t newNumCols);

  /**
   * Capacity of storage of the data structures.
   **/
  size_t cap_;

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
   * Current solution of the LP relaxation
   **/
  double *x_;

  /**
   * Current reduced costs of the variables
   **/
  double *rc_;

  /**
   * Auxiliary structure used to temporary
   * store a cut.
   **/
  OsiRowCut osrc_;

  /**
   * Lifting strategy: 0 = no lifting,
   * 1 = only one variable as wheel center,
   * 2 = a clique as wheel center
   **/
  size_t extMethod_;
};

#endif
