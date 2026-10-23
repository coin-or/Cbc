// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverHeuristics.hpp
    \brief Routines for doing heuristics.
*/

#ifndef CbcSolverHeuristics_H
#define CbcSolverHeuristics_H
#include <string>

void crunchIt(ClpSimplex *model);
class CbcHeuristicSppGrasp;
class CbcParameters;
void configureSppGraspCoverageTuning(CbcHeuristicSppGrasp &heuristic, CbcParameters &parameters,
  const std::string &stage);

/*
  On input
  doAction - 0 just fix in original and return NULL
             1 return fixed non-presolved solver
             2 as one but use presolve Inside this
             3 use presolve and fix ones with large cost
             ? do heuristics and set best solution
             ? do BAB and just set best solution
             10+ then use lastSolution and relax a few
             -2 cleanup afterwards if using 2
  On output - number fixed
*/
OsiClpSolverInterface *
fixVubs(CbcModel &model, int skipZero2,
  int &doAction,
  CoinMessageHandler * /*generalMessageHandler*/,
  const double *lastSolution, double dextra[6],
  int extra[5]);

/** 1 - add heuristics to model
        2 - do heuristics (and set cutoff and best solution)
        3 - legacy reduced-heuristics mode (skip some heuristics)
        (out model later)
    sppGraspOriginal: 1 if the original (unpreprocessed) model applies
        natively, 2 if only with artificial slacks, -1 if not, 0 unknown
    sppGraspOriginalSolution: true if the original-model call returned a
        real-feasible improving SppGrasp solution
    */
int doHeuristics(CbcModel *model, int type, CbcParameters &parameters,
  int noPrinting_, int initialPumpTune, int sppGraspOriginal = 0,
  bool sppGraspOriginalSolution = false);

#endif // CbcSolverHeuristics_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
