/* $Id: CbcSolver.hpp 1212 2009-08-21 16:19:13Z forrest $ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


/*! \file CbcSolverAnalyze.hpp
    \brief Look to see if a constraint is all-integer (variables & coeffs), or could be
  all integer.
*/

#ifndef CbcSolverAnalyze_H
#define CbcSolverAnalyze_H


int * analyze(OsiClpSolverInterface * solverMod, int & numberChanged,
		     double & increment, bool changeInt,
		     CoinMessageHandler * generalMessageHandler, bool noPrinting);
		     
#endif

