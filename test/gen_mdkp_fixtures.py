#!/usr/bin/env python3
"""
Generate Multi-Dimensional Knapsack Problem (MDKP) fixtures.

Problem:
  maximize   sum_j p[j] * x[j]
  subject to sum_j w[i][j] * x[j] <= c[i]   for i = 0..m-1
             x[j] in {0,1}                   for j = 0..n-1

Instances are solved with three independent Gurobi models to cross-validate:
  1. Standard binary MIP
  2. LP relaxation bound (to verify LP >= MIP optimal is impossible)
  3. Standard binary MIP with different branching seed

Output: mdkp_fixtures.h with embedded C arrays.
"""

import sys
import math
import random
import gurobipy as gp
from gurobipy import GRB

# --------------------------------------------------------------------------- #
# Instance specification: (m_dims, n_items, seed, alpha)                      #
# alpha controls tightness: c[i] = floor(alpha * sum_j w[i][j])              #
# --------------------------------------------------------------------------- #
SPECS = [
    (2, 20, 42,  0.50),
    (3, 25, 42,  0.50),
    (5, 30, 42,  0.50),
    (5, 40, 42,  0.50),
    (2, 20, 137, 0.50),
    (3, 25, 137, 0.50),
    (5, 30, 137, 0.50),
    (5, 40, 137, 0.50),
]


def generate_instance(m, n, seed, alpha):
    rng = random.Random(seed)
    p = [rng.randint(1, 100) for _ in range(n)]          # profits
    w = [[rng.randint(1, 30) for _ in range(n)]           # weights[dim][item]
         for _ in range(m)]
    c = [int(math.floor(alpha * sum(w[i]))) for i in range(m)]  # capacities
    return p, w, c


def solve_mdkp(m, n, p, w, c, binary=True, seed_val=0):
    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.setParam("Seed", seed_val)
    env.start()
    mdl = gp.Model(env=env)
    vtype = GRB.BINARY if binary else GRB.CONTINUOUS
    x = mdl.addVars(n, lb=0, ub=1, vtype=vtype, name="x")
    mdl.setObjective(gp.quicksum(p[j] * x[j] for j in range(n)), GRB.MAXIMIZE)
    for i in range(m):
        mdl.addConstr(gp.quicksum(w[i][j] * x[j] for j in range(n)) <= c[i],
                      name=f"cap_{i}")
    mdl.optimize()
    if mdl.Status not in (GRB.OPTIMAL,):
        raise RuntimeError(f"Model not optimal: status={mdl.Status}")
    obj = mdl.ObjVal
    mdl.close(); env.dispose()
    return obj


def arr_c(values, indent=4):
    """Format a list of ints as a C initializer list, 15 per line."""
    lines = []
    row = []
    for v in values:
        row.append(str(v))
        if len(row) == 15:
            lines.append(" " * indent + ", ".join(row) + ",")
            row = []
    if row:
        lines.append(" " * indent + ", ".join(row))
    return "\n".join(lines)


def main():
    records = []
    for (m, n, seed, alpha) in SPECS:
        tag = f"m{m:02d}_n{n:02d}_s{seed}"
        print(f"[{tag}] generating ...", flush=True)
        p, w, c = generate_instance(m, n, seed, alpha)

        # Solve 1: standard binary MIP
        obj1 = solve_mdkp(m, n, p, w, c, binary=True, seed_val=0)
        # Solve 2: binary MIP with different Gurobi seed
        obj2 = solve_mdkp(m, n, p, w, c, binary=True, seed_val=42)
        # Solve 3: LP relaxation (should be >= MIP optimal)
        lp   = solve_mdkp(m, n, p, w, c, binary=False)

        opt1, opt2 = int(round(obj1)), int(round(obj2))
        if opt1 != opt2:
            raise ValueError(f"[{tag}] cross-validation FAILED: {opt1} vs {opt2}")
        if lp < obj1 - 1e-4:
            raise ValueError(f"[{tag}] LP < MIP: {lp} < {obj1}")
        opt = opt1
        print(f"[{tag}] opt={opt}  LP_bound={lp:.2f}  cross-validated OK")
        records.append((m, n, seed, opt, p, w, c, tag))

    # ------------------------------------------------------------------ #
    # Write header                                                         #
    # ------------------------------------------------------------------ #
    lines = []
    lines.append("/* Multi-Dimensional Knapsack Problem (MDKP) test fixtures.")
    lines.append("   Profits and weights are random integers; capacities are")
    lines.append("   50% of total weight per dimension.  Optimal values cross-")
    lines.append("   validated with two independent MIP solves.")
    lines.append("   Do not edit manually; regenerate with test/gen_mdkp_fixtures.py.")
    lines.append("*/")
    lines.append("")
    lines.append("#ifndef MDKP_FIXTURES_H")
    lines.append("#define MDKP_FIXTURES_H")
    lines.append("")
    lines.append("typedef struct {")
    lines.append("  int m;           /* number of knapsack dimensions (constraints) */")
    lines.append("  int n;           /* number of items */")
    lines.append("  int opt;         /* known optimal value (maximisation) */")
    lines.append("  const int *p;    /* profits [n] */")
    lines.append("  const int *w;    /* weights [m*n], w[i*n+j] */")
    lines.append("  const int *c;    /* capacities [m] */")
    lines.append("} MdkpInstance;")
    lines.append("")

    var_lines = []
    inst_names = []

    for (m, n, seed, opt, p, w, c, tag) in records:
        comment = f"/* mdkp_{tag}: m={m} n={n} seed={seed} opt={opt} */"
        lines.append(comment)

        # profits
        pname = f"p_mdkp_{tag}"
        lines.append(f"static const int {pname}[{n}] = {{")
        lines.append(arr_c(p))
        lines.append("};")

        # weights (flattened row-major: w[i*n+j])
        wflat = [w[i][j] for i in range(m) for j in range(n)]
        wname = f"w_mdkp_{tag}"
        lines.append(f"static const int {wname}[{m*n}] = {{")
        lines.append(arr_c(wflat))
        lines.append("};")

        # capacities
        cname = f"c_mdkp_{tag}"
        lines.append(f"static const int {cname}[{m}] = {{{', '.join(map(str,c))}}};")
        lines.append("")

        inst_names.append((tag, m, n, opt, pname, wname, cname))

    # Instance array
    lines.append(f"#define NUM_MDKP_INSTANCES {len(inst_names)}")
    lines.append("static const MdkpInstance MDKP_INSTANCES[] = {")
    for (tag, m, n, opt, pname, wname, cname) in inst_names:
        lines.append(f"  /* {tag} */")
        lines.append(f"  {{ {m}, {n}, {opt}, {pname}, {wname}, {cname} }},")
    lines.append("};")
    lines.append("")
    lines.append("#endif /* MDKP_FIXTURES_H */")

    out_path = "mdkp_fixtures.h"
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nWrote {out_path} with {len(records)} instances.")


if __name__ == "__main__":
    main()
