/* $Id$ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverExpandKnapsack.hpp
    \brief Expanding possibilities of x*y, where x*y are both integers, constructing
  a knapsack constraint. Results in a tighter model.
*/

#ifndef CbcSolverExpandKnapsack_H
#define CbcSolverExpandKnapsack_H

#include "OsiClpSolverInterface.hpp"
#include "CglStored.hpp"

OsiSolverInterface *
expandKnapsack(CoinModel &model, int *whichColumn, int *knapsackStart,
  int *knapsackRow, int &numberKnapsack,
  CglStored &stored, int logLevel,
  int fixedPriority, int SOSPriority, CoinModel &tightenedModel);

void afterKnapsack(const CoinModel &coinModel2, const int *whichColumn, const int *knapsackStart,
  const int *knapsackRow, int numberKnapsack,
  const double *knapsackSolution, double *solution, int logLevel);

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
