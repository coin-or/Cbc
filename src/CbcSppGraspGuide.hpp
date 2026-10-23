// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcSppGraspGuide_H
#define CbcSppGraspGuide_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace CbcSppGraspGuide {

inline bool unchangedBounds(const std::vector< double > &snapshotLower,
  const std::vector< double > &snapshotUpper,
  const std::vector< double > &currentLower,
  const std::vector< double > &currentUpper)
{
  if (snapshotLower.size() != snapshotUpper.size()
    || snapshotLower.size() != currentLower.size()
    || snapshotLower.size() != currentUpper.size())
    return false;
  for (size_t j = 0; j < snapshotLower.size(); j++) {
    const double oldLower = snapshotLower[j];
    const double oldUpper = snapshotUpper[j];
    const double lower = currentLower[j];
    const double upper = currentUpper[j];
    if (std::isnan(oldLower) || std::isnan(oldUpper)
      || std::isnan(lower) || std::isnan(upper))
      return false;
    const auto sameBound = [](double a, double b) {
      if (a == b)
        return true;
      if (!std::isfinite(a) || !std::isfinite(b))
        return false;
      const double tolerance =
        1.0e-10 * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
      return std::fabs(a - b) <= tolerance;
    };
    if (!sameBound(oldLower, lower) || !sameBound(oldUpper, upper))
      return false;
  }
  return true;
}

inline bool preferRawSourcePoint(const std::vector< double > &rootLp,
  const std::vector< double > &pumpSourcePoint,
  const std::vector< double > &lower,
  const std::vector< double > &upper,
  const std::vector< int > &rootIntegerColumns,
  const std::vector< int > &pumpIntegerColumns,
  double pumpSourceFractionality, int pumpSourceFractionalCount,
  double integerTolerance,
  double boundTolerance, double &rootFractionality,
  int &rootFractionalCount)
{
  rootFractionality = 0.0;
  rootFractionalCount = 0;
  const size_t n = rootLp.size();
  if (!n || pumpSourcePoint.size() != n || lower.size() != n
    || upper.size() != n || rootIntegerColumns.empty()
    || rootIntegerColumns != pumpIntegerColumns
    || !std::isfinite(pumpSourceFractionality)
    || pumpSourceFractionality < 0.0 || pumpSourceFractionalCount < 0
    || !std::isfinite(integerTolerance) || integerTolerance < 0.0
    || !std::isfinite(boundTolerance) || boundTolerance < 0.0)
    return false;
  for (size_t j = 0; j < n; j++) {
    if (!std::isfinite(rootLp[j]) || !std::isfinite(pumpSourcePoint[j])
      || std::isnan(lower[j]) || std::isnan(upper[j])
      || lower[j] > upper[j]
      || rootLp[j] < lower[j] - boundTolerance
      || rootLp[j] > upper[j] + boundTolerance
      || pumpSourcePoint[j] < lower[j] - boundTolerance
      || pumpSourcePoint[j] > upper[j] + boundTolerance)
      return false;
  }
  for (size_t k = 0; k < rootIntegerColumns.size(); k++) {
    const int j = rootIntegerColumns[k];
    if (j < 0 || static_cast< size_t >(j) >= n
      || (k && rootIntegerColumns[k - 1] >= j))
      return false;
    const double value = rootLp[j];
    const double distance = std::fabs(value - std::floor(value + 0.5));
    rootFractionality += distance;
    if (distance > integerTolerance)
      rootFractionalCount++;
  }
  const double comparisonTolerance =
    1.0e-12 * std::max(1.0, rootFractionality);
  double guideFractionality = 0.0;
  int guideFractionalCount = 0;
  for (size_t k = 0; k < rootIntegerColumns.size(); k++) {
    const int j = rootIntegerColumns[k];
    const double distance =
      std::fabs(pumpSourcePoint[j] - std::floor(pumpSourcePoint[j] + 0.5));
    guideFractionality += distance;
    if (distance > integerTolerance)
      guideFractionalCount++;
  }
  const double metadataTolerance =
    1.0e-10 * std::max(1.0, pumpSourceFractionality);
  if (guideFractionalCount != pumpSourceFractionalCount
    || std::fabs(guideFractionality - pumpSourceFractionality)
      > metadataTolerance)
    return false;
  return pumpSourceFractionality
    < rootFractionality - comparisonTolerance;
}

} // namespace CbcSppGraspGuide

#endif
