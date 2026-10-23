// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcProblemClassifier.hpp"
#include "CbcParam.hpp"
#include "OsiFeatures.hpp"
#include "OsiSolverInterface.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>

const char *CbcProblemClassifier::className(ProblemClass pc)
{
  switch (pc) {
  case PureSetPartitioning:
    return "PureSetPartitioning";
  case PureSetPacking:
    return "PureSetPacking";
  case SetPartitioningPacking:
    return "SetPartitioningPacking";
  case PureSetCovering:
    return "PureSetCovering";
  case MixedSPP:
    return "MixedSPP";
  case MostlySPP:
    return "MostlySPP";
  case ZeroOneKnapsack:
    return "ZeroOneKnapsack";
  case GeneralZeroOne:
    return "GeneralZeroOne";
  case GeneralMIP:
  default:
    return "GeneralMIP";
  }
}

CbcProblemClassifier::ProblemClass CbcProblemClassifier::fromMode(
  CbcParameters::ProblemClassMode mode)
{
  switch (mode) {
  case CbcParameters::ProbClassPartitioning:
    return PureSetPartitioning;
  case CbcParameters::ProbClassPacking:
    return PureSetPacking;
  case CbcParameters::ProbClassPartPack:
    return SetPartitioningPacking;
  case CbcParameters::ProbClassCovering:
    return PureSetCovering;
  case CbcParameters::ProbClassMixedSPP:
    return MixedSPP;
  case CbcParameters::ProbClassMostlySPP:
    return MostlySPP;
  case CbcParameters::ProbClassKnapsack:
    return ZeroOneKnapsack;
  case CbcParameters::ProbClassGeneral01:
    return GeneralZeroOne;
  case CbcParameters::ProbClassGeneral:
  case CbcParameters::ProbClassAuto:
  default:
    return GeneralMIP;
  }
}

CbcProblemClassifier::ClassificationResult CbcProblemClassifier::classify(
  OsiSolverInterface *solver,
  CbcParameters::ProblemClassMode overrideMode)
{
  if (!solver) {
    ClassificationResult res;
    res.problemClass = GeneralMIP;
    res.className = className(GeneralMIP);
    res.rationale = "Null solver interface";
    return res;
  }

  std::vector< double > features(OsiFeatures::n, 0.0);
  OsiFeatures::compute(features.data(), solver);
  return classify(features.data(), overrideMode);
}

CbcProblemClassifier::ClassificationResult CbcProblemClassifier::classify(
  const double *features,
  CbcParameters::ProblemClassMode overrideMode)
{
  ClassificationResult res;
  if (!features) {
    res.problemClass = GeneralMIP;
    res.className = className(GeneralMIP);
    res.rationale = "Null features pointer";
    return res;
  }

  res.numCols = static_cast< int >(features[OFcols]);
  res.numRows = static_cast< int >(features[OFrows]);
  res.numBin = static_cast< int >(features[OFbin]);
  res.percBin = features[OFpercBin];

  double partRows = features[OFrowsPartitioning];
  double packRows = features[OFrowsPacking];
  double covRows = features[OFrowsCovering];
  double knapRows = features[OFrowsKnapsack] + features[OFrowsIntegerKnapsack];
  double sppRows = partRows + packRows + covRows;

  double totalRows = (res.numRows > 0) ? static_cast< double >(res.numRows) : 1.0;
  res.percPartitioning = (partRows / totalRows) * 100.0;
  res.percPacking = (packRows / totalRows) * 100.0;
  res.percCovering = (covRows / totalRows) * 100.0;
  res.percSPP = (sppRows / totalRows) * 100.0;
  res.percKnapsack = (knapRows / totalRows) * 100.0;

  if (overrideMode != CbcParameters::ProbClassAuto) {
    res.problemClass = fromMode(overrideMode);
    res.className = className(res.problemClass);
    res.rationale = "User override via -problemClass";
    return res;
  }

  // Pure binary check: all variables are binary
  bool isPureBinary = (res.numCols > 0 && res.numBin == res.numCols) || (res.percBin >= 99.99);

  if (isPureBinary) {
    if (res.percPartitioning >= 99.99 || (partRows == res.numRows && res.numRows > 0)) {
      res.problemClass = PureSetPartitioning;
      res.className = className(PureSetPartitioning);
      res.rationale = "100% binary variables, 100% set-partitioning constraints (sum x_j = 1)";
    } else if (res.percPacking >= 99.99 || (packRows == res.numRows && res.numRows > 0)) {
      res.problemClass = PureSetPacking;
      res.className = className(PureSetPacking);
      res.rationale = "100% binary variables, 100% set-packing constraints (sum x_j <= 1)";
    } else if (res.percCovering < 0.01 && (res.percPartitioning + res.percPacking >= 85.0 || sppRows == res.numRows)
      && (features[OFaMax] <= 1.0001)) {
      res.problemClass = SetPartitioningPacking;
      res.className = className(SetPartitioningPacking);
      res.rationale = "100% binary variables, 0-1 coefficients, partitioning and packing constraints (sum x_j = 1 or <= 1, like rmp11)";
    } else if (res.percCovering >= 99.99 || (covRows == res.numRows && res.numRows > 0)) {
      res.problemClass = PureSetCovering;
      res.className = className(PureSetCovering);
      res.rationale = "100% binary variables, 100% set-covering constraints (sum x_j >= 1)";
    } else if (res.percSPP >= 99.99 || (sppRows == res.numRows && res.numRows > 0)) {
      res.problemClass = MixedSPP;
      res.className = className(MixedSPP);
      res.rationale = "100% binary variables, 100% SPP constraints (mix of partitioning, packing, covering)";
    } else if (res.percSPP >= 85.0) {
      res.problemClass = MostlySPP;
      res.className = className(MostlySPP);
      std::ostringstream oss;
      oss << "100% binary variables, " << static_cast< int >(res.percSPP + 0.5) << "% SPP constraints";
      res.rationale = oss.str();
    } else if (res.percKnapsack >= 85.0) {
      res.problemClass = ZeroOneKnapsack;
      res.className = className(ZeroOneKnapsack);
      std::ostringstream oss;
      oss << "100% binary variables, " << static_cast< int >(res.percKnapsack + 0.5) << "% 0-1 knapsack constraints";
      res.rationale = oss.str();
    } else {
      res.problemClass = GeneralZeroOne;
      res.className = className(GeneralZeroOne);
      res.rationale = "100% binary variables, general constraint structure";
    }
  } else {
    res.problemClass = GeneralMIP;
    res.className = className(GeneralMIP);
    res.rationale = "General mixed-integer or continuous variables present";
  }

  return res;
}

const std::vector< int > &CbcProblemClassifier::autoCutParams()
{
  static const int codes[] = {
    CbcParam::CLIQUECUTS, CbcParam::FLOWCUTS, CbcParam::GMICUTS,
    CbcParam::GOMORYCUTS, CbcParam::KNAPSACKCUTS, CbcParam::LANDPCUTS,
    CbcParam::MIRCUTS, CbcParam::ODDWHEELCUTS, CbcParam::PROBINGCUTS,
    CbcParam::REDSPLITCUTS, CbcParam::REDSPLIT2CUTS, CbcParam::TWOMIRCUTS,
    CbcParam::ZEROHALFCUTS
  };
  static const std::vector< int > list(codes, codes + sizeof(codes) / sizeof(codes[0]));
  return list;
}

CbcProblemClassifier::CutTuningRecommendations CbcProblemClassifier::getCutTuning(
  const ClassificationResult &result)
{
  CutTuningRecommendations rec;

  switch (result.problemClass) {
  case PureSetPartitioning:
  case SetPartitioningPacking: {
    rec.apply = true;
    rec.summary = (result.problemClass == PureSetPartitioning)
      ? "Pure Set Partitioning: active: BKClique, OddWheel, Probing; disabled: Gomory, TwoMir, MIR2, GMI, LandP, RedSplit2, Knapsack, FlowCover, ZeroHalf"
      : "Set Partitioning & Packing: active: BKClique, OddWheel, Probing; disabled: Gomory, TwoMir, MIR2, GMI, LandP, RedSplit2, Knapsack, FlowCover, ZeroHalf";
    // Disabling generators known to produce 0 cuts on pure partitioning & packing
    rec.disableParams.push_back(CbcParam::GOMORYCUTS);
    rec.disableParams.push_back(CbcParam::TWOMIRCUTS);
    rec.disableParams.push_back(CbcParam::MIRCUTS);
    rec.disableParams.push_back(CbcParam::GMICUTS);
    rec.disableParams.push_back(CbcParam::REDSPLIT2CUTS);
    rec.disableParams.push_back(CbcParam::REDSPLITCUTS);
    rec.disableParams.push_back(CbcParam::LANDPCUTS);
    rec.disableParams.push_back(CbcParam::KNAPSACKCUTS);
    rec.disableParams.push_back(CbcParam::FLOWCUTS);
    rec.disableParams.push_back(CbcParam::ZEROHALFCUTS);

    // Active generators
    rec.enableParams.push_back(std::make_pair(CbcParam::CLIQUECUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::ODDWHEELCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::PROBINGCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    break;
  }

  case PureSetPacking: {
    rec.apply = true;
    rec.summary = "Pure Set Packing: active: BKClique, OddWheel, Gomory, TwoMir, Probing; disabled: Knapsack, FlowCover, LandP, RedSplit2, GMI";
    rec.disableParams.push_back(CbcParam::KNAPSACKCUTS);
    rec.disableParams.push_back(CbcParam::FLOWCUTS);
    rec.disableParams.push_back(CbcParam::LANDPCUTS);
    rec.disableParams.push_back(CbcParam::REDSPLIT2CUTS);
    rec.disableParams.push_back(CbcParam::REDSPLITCUTS);
    rec.disableParams.push_back(CbcParam::GMICUTS);

    rec.enableParams.push_back(std::make_pair(CbcParam::CLIQUECUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::ODDWHEELCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::GOMORYCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::TWOMIRCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::PROBINGCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    break;
  }

  case PureSetCovering: {
    rec.apply = true;
    rec.summary = "Pure Set Covering: active: ZeroHalf, Gomory, TwoMir, Probing; disabled: Clique, OddWheel, Knapsack, FlowCover, LandP, RedSplit2";
    rec.disableParams.push_back(CbcParam::CLIQUECUTS);
    rec.disableParams.push_back(CbcParam::ODDWHEELCUTS);
    rec.disableParams.push_back(CbcParam::KNAPSACKCUTS);
    rec.disableParams.push_back(CbcParam::FLOWCUTS);
    rec.disableParams.push_back(CbcParam::LANDPCUTS);
    rec.disableParams.push_back(CbcParam::REDSPLIT2CUTS);
    rec.disableParams.push_back(CbcParam::REDSPLITCUTS);
    rec.disableParams.push_back(CbcParam::GMICUTS);

    rec.enableParams.push_back(std::make_pair(CbcParam::ZEROHALFCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::GOMORYCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::TWOMIRCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::PROBINGCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    break;
  }

  case MixedSPP:
  case MostlySPP: {
    rec.apply = true;
    rec.summary = "Mixed/Mostly SPP: active: ZeroHalf, BKClique, OddWheel, Gomory, TwoMir, MIR2, Probing; disabled: Knapsack, FlowCover, LandP, RedSplit2";
    rec.disableParams.push_back(CbcParam::KNAPSACKCUTS);
    rec.disableParams.push_back(CbcParam::FLOWCUTS);
    rec.disableParams.push_back(CbcParam::LANDPCUTS);
    rec.disableParams.push_back(CbcParam::REDSPLIT2CUTS);
    rec.disableParams.push_back(CbcParam::REDSPLITCUTS);

    rec.enableParams.push_back(std::make_pair(CbcParam::ZEROHALFCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::CLIQUECUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::ODDWHEELCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::GOMORYCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::TWOMIRCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::MIRCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::PROBINGCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    break;
  }

  case ZeroOneKnapsack: {
    rec.apply = true;
    rec.summary = "0-1 Knapsack: active: Knapsack, Gomory, TwoMir, Clique, Probing";
    rec.enableParams.push_back(std::make_pair(CbcParam::KNAPSACKCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::GOMORYCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::TWOMIRCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::CLIQUECUTS, static_cast< int >(CbcParameters::CGIfMove)));
    rec.enableParams.push_back(std::make_pair(CbcParam::PROBINGCUTS, static_cast< int >(CbcParameters::CGIfMove)));
    break;
  }

  case GeneralZeroOne:
  case GeneralMIP:
  default:
    rec.apply = false;
    break;
  }

  return rec;
}
