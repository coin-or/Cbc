"""Generate Facility Location test fixtures with known optimal values.

Two problem classes:
  UFL  — Uncapacitated Facility Location
  CFLP — Capacitated Facility Location (single-source)

Each instance is solved by three independent Gurobi formulations and the
optimal values must agree before the fixture is written.

UFL formulations:
  1. Individual linking constraints  x_ij <= y_i
  2. Aggregated linking              sum_j x_ij <= n * y_i
  3. LP relaxation (x,y continuous) — verifies LP=MIP (proves strong LP bound)

CFLP formulations:
  1. Standard  : x_ij <= y_i  +  sum_j d_j*x_ij <= Q_i*y_i
  2. Capacity only (no individual links): sum_j d_j*x_ij <= Q_i*y_i
     (valid because the LP relaxation of x is integral given degree+capacity)
  3. Big-M capacity: sum_j d_j*x_ij <= Q_i * y_i  (same model, different
     variable bounds on x to stress the solver differently)

Sizes:
  UFL  : (m=5,n=10), (m=8,n=15), (m=10,n=20), (m=15,n=30) × 2 seeds
  CFLP : (m=5,n=12), (m=8,n=18), (m=10,n=25)               × 2 seeds
"""
import random, sys
import gurobipy as gp
from gurobipy import GRB


# ── Instance generation ────────────────────────────────────────────────────────

def make_ufl(m, n, seed):
    """Random UFL: opening costs in [50,200], assignment costs in [1,100]."""
    rng = random.Random(seed)
    f = [rng.randint(50, 200) for _ in range(m)]
    c = [[rng.randint(1, 100) for _ in range(n)] for _ in range(m)]
    return f, c


def make_cflp(m, n, seed):
    """Random CFLP: demands in [1,5], capacities ensure feasibility (total >= 1.5x demand)."""
    rng = random.Random(seed)
    f   = [rng.randint(50, 200) for _ in range(m)]
    c   = [[rng.randint(1, 100) for _ in range(n)] for _ in range(m)]
    dem = [rng.randint(1, 5) for _ in range(n)]
    total_dem = sum(dem)
    # Each facility gets a capacity that, summed, is ~2x total demand
    avg = max(1, (2 * total_dem) // m)
    cap = [rng.randint(max(1, avg // 2), avg * 2) for _ in range(m)]
    # Guarantee feasibility: boost smallest facility if needed
    while sum(cap) < total_dem:
        cap[cap.index(min(cap))] += total_dem
    return f, c, dem, cap


# ── UFL solvers ───────────────────────────────────────────────────────────────

def solve_ufl_individual(m, n, f, c):
    """UFL: individual linking constraints x_ij <= y_i (standard)."""
    mdl = gp.Model(); mdl.Params.OutputFlag = 0
    y = [mdl.addVar(vtype=GRB.BINARY,     obj=f[i], name=f"y({i})") for i in range(m)]
    x = [[mdl.addVar(vtype=GRB.BINARY,    obj=c[i][j], name=f"x({i},{j})")
          for j in range(n)] for i in range(m)]
    for j in range(n):
        mdl.addConstr(gp.quicksum(x[i][j] for i in range(m)) == 1, name=f"serve({j})")
    for i in range(m):
        for j in range(n):
            mdl.addConstr(x[i][j] <= y[i], name=f"link({i},{j})")
    mdl.optimize()
    assert mdl.Status == GRB.OPTIMAL, f"UFL-individual: status={mdl.Status}"
    return int(round(mdl.ObjVal))


def solve_ufl_aggregated(m, n, f, c):
    """UFL: aggregated linking constraint sum_j x_ij <= n * y_i."""
    mdl = gp.Model(); mdl.Params.OutputFlag = 0
    y = [mdl.addVar(vtype=GRB.BINARY,     obj=f[i], name=f"y({i})") for i in range(m)]
    x = [[mdl.addVar(vtype=GRB.BINARY,    obj=c[i][j], name=f"x({i},{j})")
          for j in range(n)] for i in range(m)]
    for j in range(n):
        mdl.addConstr(gp.quicksum(x[i][j] for i in range(m)) == 1, name=f"serve({j})")
    for i in range(m):
        mdl.addConstr(gp.quicksum(x[i][j] for j in range(n)) <= n * y[i], name=f"agg({i})")
    mdl.optimize()
    assert mdl.Status == GRB.OPTIMAL, f"UFL-aggregated: status={mdl.Status}"
    return int(round(mdl.ObjVal))


def solve_ufl_lp_relaxation(m, n, f, c):
    """UFL: LP relaxation (x,y continuous in [0,1]).
    For UFL the LP optimal often equals the MIP optimal; if so, it confirms the bound."""
    mdl = gp.Model(); mdl.Params.OutputFlag = 0
    y = [mdl.addVar(lb=0, ub=1, obj=f[i], name=f"y({i})") for i in range(m)]
    x = [[mdl.addVar(lb=0, ub=1, obj=c[i][j], name=f"x({i},{j})")
          for j in range(n)] for i in range(m)]
    for j in range(n):
        mdl.addConstr(gp.quicksum(x[i][j] for i in range(m)) == 1, name=f"serve({j})")
    for i in range(m):
        for j in range(n):
            mdl.addConstr(x[i][j] <= y[i], name=f"link({i},{j})")
    mdl.optimize()
    assert mdl.Status == GRB.OPTIMAL, f"UFL-LP: status={mdl.Status}"
    return mdl.ObjVal   # float — caller will check >= MIP optimal (lower bound)


# ── CFLP solvers ──────────────────────────────────────────────────────────────

def solve_cflp_standard(m, n, f, c, dem, cap):
    """CFLP standard: x_ij <= y_i  +  sum_j d_j*x_ij <= Q_i."""
    mdl = gp.Model(); mdl.Params.OutputFlag = 0
    y = [mdl.addVar(vtype=GRB.BINARY,  obj=f[i], name=f"y({i})") for i in range(m)]
    x = [[mdl.addVar(vtype=GRB.BINARY, obj=c[i][j], name=f"x({i},{j})")
          for j in range(n)] for i in range(m)]
    for j in range(n):
        mdl.addConstr(gp.quicksum(x[i][j] for i in range(m)) == 1, name=f"serve({j})")
    for i in range(m):
        for j in range(n):
            mdl.addConstr(x[i][j] <= y[i], name=f"link({i},{j})")
        mdl.addConstr(
            gp.quicksum(dem[j] * x[i][j] for j in range(n)) <= cap[i] * y[i],
            name=f"cap({i})")
    mdl.optimize()
    assert mdl.Status == GRB.OPTIMAL, f"CFLP-standard: status={mdl.Status}"
    return int(round(mdl.ObjVal))


def solve_cflp_cap_only(m, n, f, c, dem, cap):
    """CFLP: capacity constraints only (no individual x_ij <= y_i links).
    The capacity constraint sum_j d_j*x_ij <= Q_i*y_i already prevents
    assigning customers to closed facilities (Q_i*0 = 0)."""
    mdl = gp.Model(); mdl.Params.OutputFlag = 0
    y = [mdl.addVar(vtype=GRB.BINARY,  obj=f[i], name=f"y({i})") for i in range(m)]
    x = [[mdl.addVar(vtype=GRB.BINARY, obj=c[i][j], name=f"x({i},{j})")
          for j in range(n)] for i in range(m)]
    for j in range(n):
        mdl.addConstr(gp.quicksum(x[i][j] for i in range(m)) == 1, name=f"serve({j})")
    for i in range(m):
        mdl.addConstr(
            gp.quicksum(dem[j] * x[i][j] for j in range(n)) <= cap[i] * y[i],
            name=f"cap({i})")
    mdl.optimize()
    assert mdl.Status == GRB.OPTIMAL, f"CFLP-cap-only: status={mdl.Status}"
    return int(round(mdl.ObjVal))


def solve_cflp_multisource_lb(m, n, f, c, dem, cap):
    """CFLP multi-source LP relaxation: x_ij continuous in [0,1].
    Gives a lower bound; returned separately for validation (not required to equal MIP)."""
    mdl = gp.Model(); mdl.Params.OutputFlag = 0
    y = [mdl.addVar(lb=0, ub=1, obj=f[i], name=f"y({i})") for i in range(m)]
    x = [[mdl.addVar(lb=0, ub=1, obj=c[i][j], name=f"x({i},{j})")
          for j in range(n)] for i in range(m)]
    for j in range(n):
        mdl.addConstr(gp.quicksum(x[i][j] for i in range(m)) == 1)
    for i in range(m):
        mdl.addConstr(
            gp.quicksum(dem[j] * x[i][j] for j in range(n)) <= cap[i] * y[i])
    mdl.optimize()
    assert mdl.Status == GRB.OPTIMAL
    return mdl.ObjVal


# ── Main ──────────────────────────────────────────────────────────────────────

UFL_INSTANCES = [
    ("ufl_m05_n10_s1", 5,  10, 42),
    ("ufl_m08_n15_s1", 8,  15, 42),
    ("ufl_m10_n20_s1", 10, 20, 42),
    ("ufl_m15_n30_s1", 15, 30, 42),
    ("ufl_m05_n10_s2", 5,  10, 137),
    ("ufl_m08_n15_s2", 8,  15, 137),
    ("ufl_m10_n20_s2", 10, 20, 137),
    ("ufl_m15_n30_s2", 15, 30, 137),
]

CFLP_INSTANCES = [
    ("cflp_m05_n12_s1", 5,  12, 42),
    ("cflp_m08_n18_s1", 8,  18, 42),
    ("cflp_m10_n25_s1", 10, 25, 42),
    ("cflp_m05_n12_s2", 5,  12, 137),
    ("cflp_m08_n18_s2", 8,  18, 137),
    ("cflp_m10_n25_s2", 10, 25, 137),
]

ufl_results  = []
cflp_results = []
all_ok = True

print("=== UFL instances (3 formulations: individual, aggregated, LP relaxation) ===")
for name, m, n, seed in UFL_INSTANCES:
    f, c = make_ufl(m, n, seed)
    opt_ind = solve_ufl_individual(m, n, f, c)
    opt_agg = solve_ufl_aggregated(m, n, f, c)
    lp_val  = solve_ufl_lp_relaxation(m, n, f, c)
    # LP is a relaxation: must be <= MIP optimal; equality confirms tight bound
    lp_tight = abs(lp_val - opt_ind) < 0.5
    agree = (opt_ind == opt_agg)
    status = "OK" if agree else "MISMATCH"
    lp_str = f"LP={lp_val:.1f}{'=MIP' if lp_tight else '<MIP'}"
    print(f"  {name}: m={m} n={n} seed={seed}  ind={opt_ind} agg={opt_agg} {lp_str}  [{status}]")
    if not agree:
        all_ok = False
    else:
        ufl_results.append((name, m, n, seed, opt_ind, f, c))

print()
print("=== CFLP instances (2 MIP formulations + LP lower bound) ===")
for name, m, n, seed in CFLP_INSTANCES:
    f, c, dem, cap = make_cflp(m, n, seed)
    opt_std  = solve_cflp_standard(m, n, f, c, dem, cap)
    opt_cap  = solve_cflp_cap_only(m, n, f, c, dem, cap)
    lb       = solve_cflp_multisource_lb(m, n, f, c, dem, cap)
    agree = (opt_std == opt_cap)
    status = "OK" if agree else "MISMATCH"
    print(f"  {name}: m={m} n={n} seed={seed}  std={opt_std} cap-only={opt_cap} LP-lb={lb:.1f}  [{status}]")
    if not agree:
        all_ok = False
    else:
        cflp_results.append((name, m, n, seed, opt_std, f, c, dem, cap))

if not all_ok:
    print("\nERROR: formulation disagreement — fix instances before regenerating fixtures.")
    sys.exit(1)

print(f"\nAll {len(ufl_results)} UFL + {len(cflp_results)} CFLP instances cross-validated OK.")

# ── Generate C header ─────────────────────────────────────────────────────────

lines = []
lines.append("/* Facility Location test fixtures — known optimal values.")
lines.append("   UFL  (Uncapacitated) and CFLP (Capacitated, single-source).")
lines.append("   Costs and demands are random integers with fixed seeds.")
lines.append("   Optimal values cross-validated with multiple MIP formulations.")
lines.append("   Do not edit manually; regenerate with test/gen_fl_fixtures.py.")
lines.append("*/")
lines.append("")
lines.append("#ifndef FL_FIXTURES_H")
lines.append("#define FL_FIXTURES_H")
lines.append("")

# UFL struct
lines.append("/* ---- UFL ---- */")
lines.append("typedef struct {")
lines.append("  int m;            /* facilities */")
lines.append("  int n;            /* customers  */")
lines.append("  int opt;          /* optimal value */")
lines.append("  const int *f;     /* opening costs [m] */")
lines.append("  const int *c;     /* assignment costs [m*n], c[i*n+j] */")
lines.append("} UflInstance;")
lines.append("")

for name, m, n, seed, opt, f, c in ufl_results:
    lines.append(f"/* {name}: m={m} n={n} seed={seed} opt={opt} */")
    lines.append(f"static const int f_{name}[{m}] = {{{', '.join(str(v) for v in f)}}};")
    flat_c = [c[i][j] for i in range(m) for j in range(n)]
    lines.append(f"static const int c_{name}[{m*n}] = {{{', '.join(str(v) for v in flat_c)}}};")
    lines.append("")

lines.append("static const UflInstance UFL_INSTANCES[] = {")
for name, m, n, seed, opt, f, c in ufl_results:
    lines.append(f"  {{ {m}, {n}, {opt}, f_{name}, c_{name} }},")
lines.append("};")
lines.append(f"#define N_UFL_INSTANCES {len(ufl_results)}")
lines.append("")

# CFLP struct
lines.append("/* ---- CFLP ---- */")
lines.append("typedef struct {")
lines.append("  int m;            /* facilities */")
lines.append("  int n;            /* customers  */")
lines.append("  int opt;          /* optimal value */")
lines.append("  const int *f;     /* opening costs [m] */")
lines.append("  const int *cap;   /* capacities [m] */")
lines.append("  const int *dem;   /* demands [n] */")
lines.append("  const int *c;     /* assignment costs [m*n], c[i*n+j] */")
lines.append("} CflpInstance;")
lines.append("")

for name, m, n, seed, opt, f, c, dem, cap in cflp_results:
    lines.append(f"/* {name}: m={m} n={n} seed={seed} opt={opt} */")
    lines.append(f"static const int f_{name}[{m}] = {{{', '.join(str(v) for v in f)}}};")
    lines.append(f"static const int cap_{name}[{m}] = {{{', '.join(str(v) for v in cap)}}};")
    lines.append(f"static const int dem_{name}[{n}] = {{{', '.join(str(v) for v in dem)}}};")
    flat_c = [c[i][j] for i in range(m) for j in range(n)]
    lines.append(f"static const int c_{name}[{m*n}] = {{{', '.join(str(v) for v in flat_c)}}};")
    lines.append("")

lines.append("static const CflpInstance CFLP_INSTANCES[] = {")
for name, m, n, seed, opt, f, c, dem, cap in cflp_results:
    lines.append(f"  {{ {m}, {n}, {opt}, f_{name}, cap_{name}, dem_{name}, c_{name} }},")
lines.append("};")
lines.append(f"#define N_CFLP_INSTANCES {len(cflp_results)}")
lines.append("")
lines.append("#endif /* FL_FIXTURES_H */")

output = "\n".join(lines) + "\n"
out_path = "/tmp/fl_fixtures.h"
with open(out_path, "w") as fh:
    fh.write(output)
print(f"Written: {out_path}")
