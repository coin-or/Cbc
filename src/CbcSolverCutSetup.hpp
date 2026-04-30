// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverCutSetup.hpp
    \brief Routine for installing cut generators on a CbcModel.
*/

#ifndef CbcSolverCutSetup_H
#define CbcSolverCutSetup_H

#include "CoinBronKerbosch.hpp"

class CbcModel;
class CbcParameters;

/// Register all cut generators on babModel based on parameter settings,
/// then apply per-generator tuning (switches, accuracy, timing, cutDepth).
void installCutGenerators(
  CbcModel &babModel,
  CbcParameters &parameters,
  int complicatedInteger,
  bool dominatedCuts,
  bool miplib,
  const std::string &cgraphMode,
  int oldCliqueMode,
  int maxCallsBK,
  int bkClqExtMethod,
  CoinBronKerbosch::PivotingStrategy bkPivotingStrategy,
  int oddWExtMethod,
  int mixedRoundStrategy);

#endif // CbcSolverCutSetup_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
