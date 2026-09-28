// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverCutSetup.hpp
    \brief Routine for installing cut generators on a CbcModel.
*/

#ifndef CbcSolverCutSetup_H
#define CbcSolverCutSetup_H

#include <map>
#include <string>

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
  const std::string &cgraphMode,
  int oldCliqueMode,
  int maxCallsBK,
  int bkClqExtMethod,
  CoinBronKerbosch::PivotingStrategy bkPivotingStrategy,
  int oddWExtMethod,
  int mixedRoundStrategy,
  std::string *switchOffChoice = NULL);

/** A parsed -cutSwitchOff specification.

    The value given to CbcCutGenerator::setSwitchOffIfLessThan() for each
    generator: `auto' keeps the built-in value, an integer replaces it.
    Per-generator entries win over the global one whatever their order. */
struct CbcCutSwitchOff {
  CbcCutSwitchOff()
    : allAuto(true)
    , all(0)
  {
  }
  bool allAuto;
  int all;
  /// Lower-case generator key -> value; autoValue means built-in.
  std::map< std::string, int > byKey;
  static const int autoValue;
};

/** Parse a -cutSwitchOff specification: a comma-separated list of
    `VALUE' (every generator) and `NAME:VALUE' (one generator; `NAME=VALUE'
    is also accepted, but not on the command line) items, where
    VALUE is `auto' or an integer >= -2 and NAME is a cut option name
    without its "Cuts" suffix (e.g. twoMir, zeroHalf), case-insensitive.
    Returns false, with a message in *error, if the specification is
    malformed. */
bool parseCutSwitchOff(const std::string &spec, CbcCutSwitchOff &result,
  std::string *error);

#endif // CbcSolverCutSetup_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
