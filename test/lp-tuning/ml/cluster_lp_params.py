import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import gmean
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform, jaccard
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold

EXP_DIR      = "/home/haroldo/experiments/cbc/lp_relax_2026_05_15_noblas"
FEATURES_CSV = "/home/haroldo/inst/miplib/2017+spp/features.csv"
TIMELIMIT    = 10800.0; PENALTY = TIMELIMIT * 2.0; BASELINE = "dual_default"; N_SPLITS = 5

DROP_PARAMS = {
    "primal_PEsteep_idiot100","primal_idiot200_pertv61",
    "primal_idiot200_pertvm1483","primal_idiot200_scaling_equi",
    "primal_idiot300_pertvm1483","primal_idiot40_pertvm1483",
}

feat_df  = pd.read_csv(FEATURES_CSV).rename(columns={"Name":"instance"}).set_index("instance")
times_df = pd.read_csv(f"{EXP_DIR}/lp_avg_times_full.csv")
times_df = times_df[~times_df["param_tag"].isin(DROP_PARAMS)]
pivot    = times_df.pivot(index="instance", columns="param_tag",
                          values="avg_wall_seconds").fillna(PENALTY)

params       = list(pivot.columns)
feature_cols = list(feat_df.columns)
common       = sorted(set(feat_df.index) & set(pivot.index))
M            = pivot.loc[common].values   # (n_inst, n_params)
print(f"Instances: {len(common)}, Params: {len(params)}")

def build_linkage(threshold_pct):
    """
    Win profile: param p is 'good' on instance i if it is within threshold_pct
    of the best param on that instance (and solved, i.e. < PENALTY).
    Cluster params by Jaccard distance on these binary profiles.
    """
    best_t  = M.min(axis=1, keepdims=True)
    solved  = (M < PENALTY)
    good    = solved & (M <= best_t * (1 + threshold_pct / 100.0))  # (n_inst, n_params)
    W = good.T.astype(float)  # (n_params, n_inst)
    n = len(params)
    dist_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = jaccard(W[i], W[j])
            dist_mat[i, j] = dist_mat[j, i] = d
    np.fill_diagonal(dist_mat, 0.0)
    return linkage(squareform(dist_mat), method="ward"), good

def get_reps(Z, good, k):
    """
    For each cluster, pick the rep = param that is most often THE BEST
    (within the cluster's specialised instances), with geo-mean as tiebreaker.
    """
    labels = fcluster(Z, k, criterion="maxclust")
    reps, cluster_info = [], {}
    for c in range(1, k+1):
        midx = [i for i, l in enumerate(labels) if l == c]
        if not midx: continue
        # Instances where any cluster member is good
        spec_inst = np.where(good[:, midx].any(axis=1))[0]
        if len(spec_inst) == 0:
            spec_inst = np.arange(len(common))
        sub_M = M[np.ix_(spec_inst, midx)]  # times on specialised instances
        # Count how often each member is the fastest within the cluster on each instance
        best_in_cluster = sub_M.argmin(axis=1)         # (n_spec_inst,)
        win_counts = np.bincount(best_in_cluster, minlength=len(midx))
        # Also compute geo-mean on specialised instances (for tiebreaking & display)
        local_gm = np.expm1(np.log1p(sub_M).mean(axis=0))
        # Rep = most wins within cluster; tiebreak by local geo-mean
        best_local = sorted(range(len(midx)),
                            key=lambda j: (-win_counts[j], local_gm[j]))[0]
        rep = params[midx[best_local]]
        reps.append(rep)
        cluster_info[c] = (midx, spec_inst, win_counts, local_gm, best_local)
    return reps, cluster_info, labels

def run_kfold(reps, depth=8, n_est=150):
    sel = list(reps) if BASELINE in reps else list(reps) + [BASELINE]
    piv_sel = pivot.loc[common, sel]
    X = feat_df.loc[common][feature_cols].values.astype(np.float64)
    col_med = np.nanmedian(X, axis=0); nm = np.isnan(X)
    X[nm] = np.take(col_med, np.where(nm)[1])
    param_list = list(piv_sel.columns)
    base_idx   = param_list.index(BASELINE)
    y_labels   = piv_sel.loc[common].idxmin(axis=1).values
    le = LabelEncoder(); le.fit(y_labels); y = le.transform(y_labels)
    kf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=42)
    speedups = []
    for tr, te in kf.split(X, y):
        clf = RandomForestClassifier(n_estimators=n_est, max_depth=depth,
                                     random_state=42, n_jobs=-1)
        clf.fit(X[tr], y[tr])
        preds = le.inverse_transform(clf.predict(X[te]))
        at = [piv_sel.loc[common[i], p] if p in param_list else PENALTY
              for i, p in zip(te, preds)]
        bt = piv_sel.values[te, base_idx]
        speedups.append(float(gmean(np.clip(bt / np.maximum(at, 1e-3), 0.05, 100))))
    return np.mean(speedups), speedups

# ── Sweep threshold × k ───────────────────────────────────────────────────────
print("\n=== Competitiveness threshold × k clusters (depth=8) ===")
print(f"  baseline: 64 params d=10 → 1.298x | 24 params d=8 → 1.292x\n")

overall_best = (0, None, None, None)  # (speedup, threshold, k, reps)
for thr in [5, 10, 20]:
    Z, good = build_linkage(thr)
    print(f"  -- threshold={thr}% --")
    for k in range(4, 13):
        reps, info, _ = get_reps(Z, good, k)
        sp, folds = run_kfold(reps)
        marker = ""
        if sp > overall_best[0]:
            overall_best = (sp, thr, k, reps)
            marker = " ◀ BEST"
        print(f"    thr={thr}%  k={k:2d}  params={len(reps):2d}  "
              f"speedup={sp:.3f}x  folds={[f'{s:.2f}' for s in folds]}  {marker}")
    print()

best_sp, best_thr, best_k, best_reps = overall_best
print(f"\n=== Overall best: threshold={best_thr}%, k={best_k}, speedup={best_sp:.3f}x ===")
print(f"Reps: {best_reps}")

# ── Depth tuning on best config ───────────────────────────────────────────────
print(f"\n=== Depth tuning (threshold={best_thr}%, k={best_k}) ===")
for depth in [6, 8, 10, 12]:
    sp, folds = run_kfold(best_reps, depth=depth)
    print(f"  d={depth:2d}  speedup={sp:.3f}x  folds={[f'{s:.2f}' for s in folds]}")

# ── Print cluster contents for best config ────────────────────────────────────
Z_best, good_best = build_linkage(best_thr)
reps_best, info_best, labels_best = get_reps(Z_best, good_best, best_k)
print(f"\n=== Cluster contents (threshold={best_thr}%, k={best_k}) ===")
global_gm = np.expm1(np.log1p(M).mean(axis=0))
for c, (midx, spec_inst, win_counts, local_gm, best_local) in sorted(info_best.items()):
    rep = params[midx[best_local]]
    members = sorted(zip([params[i] for i in midx],
                         [global_gm[i] for i in midx],
                         [good_best[:, i].sum() for i in midx],
                         win_counts),
                     key=lambda x: -x[3])  # sort by wins-within-cluster
    excl = good_best[:, midx].any(axis=1) & \
           ~np.delete(good_best, midx, axis=1).any(axis=1)
    print(f"\n  Cluster {c}  ({len(midx)} params | "
          f"{len(spec_inst)} competitive instances | "
          f"{excl.sum()} exclusive) → rep: {rep}")
    print(f"    {'param':<52s} {'global-gm':>10} {'good-inst':>9} {'cluster-wins':>12}")
    for p, ggm, n_good, wc in members:
        mark = " ★" if p == rep else ""
        print(f"    {p:<52s} {ggm:10.2f}s {n_good:9d} {wc:12d}{mark}")
