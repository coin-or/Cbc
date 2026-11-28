// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcProblemClassifier_H
#define CbcProblemClassifier_H

#include "CbcConfig.h"
#include "CbcParameters.hpp"
#include <string>
#include <utility>
#include <vector>

class OsiSolverInterface;

/*! \class CbcProblemClassifier
    \brief Fast, non-destructive instance classification and algorithm tuning recommendations.

    Uses OsiFeatures (computed in O(nz) time) to categorize MIP instances into
    structural families:
      - PureSetPartitioning: 100% binary, 100% partitioning rows (sum x_j = 1)
      - PureSetPacking:      100% binary, 100% packing rows (sum x_j <= 1)
      - PureSetCovering:     100% binary, 100% covering rows (sum x_j >= 1)
      - MixedSPP:            100% binary, 100% SPP rows (partitioning + packing + covering)
      - MostlySPP:           100% binary, >= 85% SPP rows
      - ZeroOneKnapsack:     100% binary, >= 85% 0-1 knapsack rows
      - GeneralZeroOne:      100% binary, general constraint structure
      - GeneralMIP:          General mixed-integer or continuous variables present

    Provides targeted cut generator selection and filtering recommendations,
    automatically disabling ineffective/barren generators that waste root time
    while preserving or enhancing dual bounds.
*/
class CbcProblemClassifier {
public:
  enum ProblemClass {
    GeneralMIP = 0,
    PureSetPartitioning,
    PureSetPacking,
    SetPartitioningPacking,
    PureSetCovering,
    MixedSPP,
    MostlySPP,
    ZeroOneKnapsack,
    GeneralZeroOne
  };

  struct ClassificationResult {
    ProblemClass problemClass;
    std::string className;
    std::string rationale;
    double percBin;
    double percPartitioning;
    double percPacking;
    double percCovering;
    double percSPP;
    double percKnapsack;
    int numCols;
    int numRows;
    int numBin;
  };

  struct CutTuningRecommendations {
    bool apply;
    std::string summary;
    // Cut generator parameter codes to turn OFF if at default
    std::vector< int > disableParams;
    // Cut generator parameter codes to set to specific mode if at default
    std::vector< std::pair< int, int > > enableParams;
    // Score mode for cut filter (-1 if keep default)
    int filterScoreMode;
    // Whether to force filtering on all rounds
    bool filterAlways;
    // Pass and drop multipliers
    double cutPassMultiplier;
    double minDropMultiplier;

    CutTuningRecommendations()
      : apply(false)
      , filterScoreMode(-1)
      , filterAlways(false)
      , cutPassMultiplier(1.0)
      , minDropMultiplier(1.0)
    {
    }
  };

  /// Classify instance directly from solver
  static ClassificationResult classify(
    OsiSolverInterface *solver,
    CbcParameters::ProblemClassMode overrideMode = CbcParameters::ProbClassAuto);

  /// Classify from precomputed OsiFeatures vector
  static ClassificationResult classify(
    const double *features,
    CbcParameters::ProblemClassMode overrideMode = CbcParameters::ProbClassAuto);

  /// Compute cut tuning recommendations based on classification result
  static CutTuningRecommendations getCutTuning(const ClassificationResult &result);

  /// Cut generator parameters that may be `auto' (resolved by the tuner)
  static const std::vector< int > &autoCutParams();

  /// Map ProblemClass enum to readable name
  static const char *className(ProblemClass pc);

  /// Map CbcParameters::ProblemClassMode to ProblemClass
  static ProblemClass fromMode(CbcParameters::ProblemClassMode mode);
};

#endif
