/* $Id: CbcSolverExpandKnapsack.hpp 1212 2009-08-21 16:19:13Z forrest $ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


/*! \file CbcSolverExpandKnapsack.hpp
    \brief Expanding possibilities of x*y, where x*y are both integers, constructing
  a knapsack constraint. Results in a tighter model.
*/


#ifndef CbcSolverExpandKnapsack_H
#define CbcSolverExpandKnapsack_H

OsiSolverInterface *
expandKnapsack(CoinModel & model, int * whichColumn, int * knapsackStart,
               int * knapsackRow, int &numberKnapsack,
               CglStored & stored, int logLevel,
               int fixedPriority, int SOSPriority, CoinModel & tightenedModel);

#endif

