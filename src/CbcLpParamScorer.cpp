// Copyright (C) 2026, COIN-OR Foundation
// This file is distributed under the Eclipse Public License.

#include "CbcLpParamScorer.hpp"
#include "OsiFeatures.hpp"
#include <cmath>

namespace {
double feature(const double *features, OsiFeature index, double missing)
{
  return std::isnan(features[index]) ? missing : features[index];
}
}

const char *cbcRecommendLpParam(const double *features)
{
  // First match wins. Each rule retains its discovery-time NaN imputations.
  // See test/lp-tuning/doc/sequential-auto-policy.md for provenance.
  if (feature(features, OFnzPercRowsCovering, 0.0) > 50.7092105993763
    && feature(features, OFcolsLeast4Nz, 3271.5) > 1056.5
    && feature(features, OFpercRowsLeast4096Nz, 0.0) <= 0.0445202274795385)
    return "primal_idiot60";

  if (feature(features, OFnzPercRowsHubImplication, 0.0) > 55.30788084035105
    && feature(features, OFpercBin, 99.52230035549104) > 99.04976740485506
    && feature(features, OFcols, 10006.0) > 810.0)
    return "primal_idiot50";

  if (feature(features, OFnzpercRowsPartPacking, 6.06743448506093) > 99.82268606961816
    && feature(features, OFrowsLeast8Nz, 639.0) > 1135.5)
    return "primal_idiot60";

  if (feature(features, OFobjStdDev, 2.2684619995820503) <= 0.001552965637549735
    && feature(features, OFnzrowsPartPacking, 2170.5) > 2026.0)
    return "primal_idiot50_pertv50";

  return "dual_pertv72";
}
