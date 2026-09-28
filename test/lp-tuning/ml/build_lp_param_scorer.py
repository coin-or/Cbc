#!/usr/bin/env python3
"""
build_lp_param_scorer.py — Train the LP parameter recommendation model on all
available instances and emit the C++ scorer source files.

Outputs:
  Cbc/src/CbcLpParamScorer.hpp  — public API header
  Cbc/src/CbcLpParamScorer.cpp  — generated scorer implementation

The scorer is a Random Forest classifier trained on OsiFeatures (207 features)
to predict the best LP parameter configuration per instance.  At inference time
the generated C++ function `cbcRecommendLpParam(const double* features)` returns
a const string with the recommended parameter tag (e.g. "dual_pesteep_psineg1").

Model selection history (k-fold cross-validation, speedup vs dual_default):
  Round 1 — 24 params, 380 instances, timelimit=10800s:
    RF  n=100, depth=8  → 1.316× (22K lines)
    RF  n=150, depth=8  → 1.334× (34K lines)
    RF  n=200, depth=8  → 1.311× (45K lines)
    LGBM lr=0.05, leaves=31, n=200 → 1.322× (283K lines — too large)
    XGBoost 3.x not supported by m2cgen 0.10.

  Round 2 — 64 params (70 extended − 6 dropped), 380 instances, timelimit=10800s:
    RF  n=150, depth=8   → 1.186× (31K lines)  ← worse due to 64-class problem
    RF  n=150, depth=10  → 1.298× (47K lines)
    RF  n=200, depth=10  → 1.289×
    RF  n=300, depth=8   → 1.202×
    24-param set, depth=8  → 1.292× (reference on same data/timelimit)
    64-param set, depth=10 → 1.298× ← marginally better on same data

  Round 3 — 12 cluster reps (Jaccard win-profile clustering, threshold=5%, k=12),
             380 instances, timelimit=10800s:
    RF  n=150, depth=8  → 1.362× (40K lines)  ← CHOSEN — best result
    Cluster representatives:
      dual_pertv72, dual_pesteep_psi1_pertv61, dual_pesteep_scaling_off,
      primal_idiot30, dual_pesteep_psineg1_pertv61, primal_sprint,
      dual_pertv58, primal_idiot10, dual_pesteep_pertv58,
      primal_idiot50, primal_idiot500, primal_idiot60
    Speedup per fold: 1.32, 1.31, 1.42, 1.43, 1.32 — no regressions

  Why clustering helps: Jaccard distance on binary "win or competitive" profiles
  (threshold=5%) groups params that specialise on the same instances.  Picking
  one rep per cluster gives 12 classes instead of 64 — fewer classes means each
  class has more training examples, allowing shallower trees (depth=8 sufficient)
  and better generalisation.

  Dropped params (never win, competitive on ≤5 instances):
    primal_PEsteep_idiot100, primal_idiot200_pertv61,
    primal_idiot200_pertvm1483, primal_idiot200_scaling_equi,
    primal_idiot300_pertvm1483, primal_idiot40_pertvm1483

RandomForest n=150, depth=8 with 12 cluster reps is the current best model.
To recluster / retrain, see lp_tune_scripts/cluster_lp_params.py.
"""

import os
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
import m2cgen as m2c

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
REPO_DIR     = os.path.dirname(SCRIPT_DIR)
EXP_DIR      = os.path.expanduser(
    "~/experiments/cbc/lp_relax_2026_05_15_noblas")
FEATURES_CSV = "/home/haroldo/inst/miplib/2017+spp/features.csv"
AVG_TIMES    = os.path.join(EXP_DIR, "lp_avg_times_full.csv")
OUT_HPP      = os.path.join(REPO_DIR, "Cbc", "src", "CbcLpParamScorer.hpp")
OUT_CPP      = os.path.join(REPO_DIR, "Cbc", "src", "CbcLpParamScorer.cpp")

# ---------------------------------------------------------------------------
# 12 cluster representatives (Jaccard win-profile, threshold=5%, k=12)
# See lp_tune_scripts/cluster_lp_params.py for how these were selected.
# ---------------------------------------------------------------------------
REPS = [
    'dual_pertv72', 'dual_pesteep_psi1_pertv61', 'dual_pesteep_scaling_off',
    'primal_idiot30', 'dual_pesteep_psineg1_pertv61', 'primal_sprint',
    'dual_pertv58', 'primal_idiot10', 'dual_pesteep_pertv58',
    'primal_idiot50', 'primal_idiot500', 'primal_idiot60',
]

# ---------------------------------------------------------------------------
# Best hyperparameters (k-fold validated, round 3)
# RF n=150, depth=8 with 12 cluster reps: 1.362× speedup, ~40K generated C lines
# ---------------------------------------------------------------------------
RF_PARAMS = dict(
    n_estimators = 150,
    max_depth    = 8,
    random_state = 42,
    n_jobs       = -1,
)
TIMELIMIT  = 10800.0
PENALTY    = 2.0

# ---------------------------------------------------------------------------
# Load and prepare data
# ---------------------------------------------------------------------------
print("Loading features …")
feat_df = pd.read_csv(FEATURES_CSV)
feat_df = feat_df.rename(columns={"Name": "instance"})

print("Loading avg times …")
times_df = pd.read_csv(AVG_TIMES)

# Pivot to wide: rows=instance, columns=param_tag, values=avg_wall_seconds
pivot = times_df.pivot(index="instance", columns="param_tag",
                       values="avg_wall_seconds")

# Filter to cluster representatives only (adds dual_default for baseline if missing)
BASELINE = "dual_default"
sel = list(REPS) if BASELINE in REPS else REPS + [BASELINE]
missing = [p for p in sel if p not in pivot.columns]
if missing:
    print(f"WARNING: params not found in avg_times: {missing}")
    sel = [p for p in sel if p in pivot.columns]
pivot = pivot[sel].fillna(TIMELIMIT * PENALTY)
params = pivot.columns.tolist()
print(f"Parameters ({len(params)}): {params}")

# Align features to instances present in both tables
common = sorted(set(feat_df["instance"]) & set(pivot.index))
print(f"Common instances: {len(common)}")

feat_df = feat_df.set_index("instance").loc[common]
pivot   = pivot.loc[common]

feature_cols = [c for c in feat_df.columns]
print(f"Features: {len(feature_cols)}")

X = feat_df[feature_cols].values.astype(np.float64)
# Label: param with minimum average penalised time for each instance
y_labels = pivot.idxmin(axis=1).values

# Global label encoder (all instances, so every class is seen)
le = LabelEncoder()
le.fit(y_labels)
y = le.transform(y_labels)
classes = le.classes_.tolist()
n_classes = len(classes)
print(f"Classes ({n_classes}): {classes}")

# ---------------------------------------------------------------------------
# NaN imputation: fill with column median
# ---------------------------------------------------------------------------
col_medians = np.nanmedian(X, axis=0)
nan_mask    = np.isnan(X)
X[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])

# ---------------------------------------------------------------------------
# Train on ALL instances
# ---------------------------------------------------------------------------
print("Training RandomForest on all instances …")
clf = RandomForestClassifier(**RF_PARAMS)
clf.fit(X, y)
print("Training complete.")

# Quick sanity: training accuracy
train_pred    = clf.predict(X)
train_correct = (train_pred == y).mean()
print(f"Training accuracy: {train_correct:.3f}")

# ---------------------------------------------------------------------------
# Export to C via m2cgen
# ---------------------------------------------------------------------------
print("Exporting model to C via m2cgen …")
c_code = m2c.export_to_c(clf, function_name="cbcLpParamScore_impl")
print(f"Generated C code: {len(c_code):,} chars, "
      f"{c_code.count(chr(10)):,} lines")

# ---------------------------------------------------------------------------
# Feature names (must match OsiFeatures::name(i) order)
# ---------------------------------------------------------------------------
feat_names_literal = "{\n" + ",\n".join(
    f'    "{name}"' for name in feature_cols) + "\n}"

# ---------------------------------------------------------------------------
# Class names as C string literals
# ---------------------------------------------------------------------------
class_names_literal = "{\n" + ",\n".join(
    f'    "{name}"' for name in classes) + "\n}"

# ---------------------------------------------------------------------------
# Write CbcLpParamScorer.hpp — public API only
# ---------------------------------------------------------------------------
hpp_content = """\
// CbcLpParamScorer.hpp — LP parameter recommendation for MIPster
//
// Auto-generated by lp_tune_scripts/build_lp_param_scorer.py
// DO NOT EDIT MANUALLY — re-run the script to regenerate.
//
// Model: RandomForest classifier (n_estimators=150, max_depth=8)
// Training data: {n_inst} instances x {n_feat} OsiFeatures
// Classes ({n_cls}): {cls_preview}
#pragma once

/// Returns the recommended LP parameter tag for the given feature vector.
/// @param features  Array of OsiFeatures::OFCount doubles (see OsiFeatures.hpp).
/// @return  A const C string with the recommended parameter tag, e.g.
///          "dual_pesteep_psineg1".  The pointer is valid for the lifetime
///          of the program (points into a static table).
const char *cbcRecommendLpParam(const double *features);
""".format(
    n_inst=len(common),
    n_feat=len(feature_cols),
    n_cls=n_classes,
    cls_preview=", ".join(classes[:6]) + (", …" if n_classes > 6 else ""),
)

print(f"Writing {OUT_HPP} …")
with open(OUT_HPP, "w") as fh:
    fh.write(hpp_content)

# ---------------------------------------------------------------------------
# Write CbcLpParamScorer.cpp — generated scorer + public wrapper
# ---------------------------------------------------------------------------
cpp_preamble = """\
// CbcLpParamScorer.cpp — LP parameter recommendation for MIPster
//
// Auto-generated by lp_tune_scripts/build_lp_param_scorer.py
// DO NOT EDIT MANUALLY — re-run the script to regenerate.
//
// Model: RandomForest classifier (n_estimators=150, max_depth=8)
// Training set: {n_inst} instances, {n_feat} features, {n_cls} classes
#include "CbcLpParamScorer.hpp"
#include <cstring>
#include <cmath>
#include <cfloat>

// ---------------------------------------------------------------------------
// Feature metadata (informational — order matches OsiFeatures::name(i))
// ---------------------------------------------------------------------------
static const int kNFeatures = {n_feat};
static const int kNClasses  = {n_cls};

static const char *const kFeatureNames[{n_feat}] = {feat_names};

static const char *const kClassNames[{n_cls}] = {class_names};

// ---------------------------------------------------------------------------
// m2cgen-generated scorer (XGBoost per-class scores)
// Signature: cbcLpParamScore_impl(double* input, double* output)
//   input  — feature vector of length kNFeatures
//   output — score vector of length kNClasses; argmax gives predicted class
// ---------------------------------------------------------------------------
""".format(
    n_inst=len(common),
    n_feat=len(feature_cols),
    n_cls=n_classes,
    feat_names=feat_names_literal,
    class_names=class_names_literal,
)

cpp_suffix = """
// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
const char *cbcRecommendLpParam(const double *features)
{
  double scores[kNClasses];
  cbcLpParamScore_impl(const_cast<double *>(features), scores);
  int best = 0;
  for (int i = 1; i < kNClasses; ++i)
    if (scores[i] > scores[best])
      best = i;
  return kClassNames[best];
}
"""

print(f"Writing {OUT_CPP} …")
with open(OUT_CPP, "w") as fh:
    fh.write(cpp_preamble)
    fh.write(c_code)
    fh.write(cpp_suffix)

cpp_lines = cpp_preamble.count("\n") + c_code.count("\n") + cpp_suffix.count("\n")
print(f"Wrote {OUT_CPP}: {cpp_lines:,} lines")
print("Done.")
