// Copyright (C) 2026, COIN-OR Foundation
// This file is distributed under the Eclipse Public License.
#pragma once

/// Select root-only LP settings using four ordered exceptions to dual_pertv72.
/// @param features Non-null array of OFCount doubles from OsiFeatures::compute,
///                 after root strengthening and before Clp LP presolve.
/// @return A static parameter tag; no concurrent racing is requested here.
const char *cbcRecommendLpParam(const double *features);
