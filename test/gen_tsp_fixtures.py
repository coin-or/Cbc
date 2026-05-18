"""Generate random Euclidean TSP fixtures and solve optimally with Gurobi.

Three formulations are used and cross-validated per instance:
  1. DFJ-lazy   : undirected edges, subtour elimination via lazy callbacks (MIPSOL)
  2. DFJ-cuts   : undirected edges, subtour elimination via user cuts on LP relaxation
  3. MTZ         : directed arcs, compact Miller-Tucker-Zemlin formulation (no callbacks)
All three must agree on the optimal value; any disagreement aborts generation.
"""
import random, math, sys
from itertools import combinations
import gurobipy as gp
from gurobipy import GRB


def make_instance(n, seed):
    rng = random.Random(seed)
    coords = [(rng.randint(0, 200), rng.randint(0, 200)) for _ in range(n)]
    dist = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                dx = coords[i][0] - coords[j][0]
                dy = coords[i][1] - coords[j][1]
                dist[i][j] = round(math.sqrt(dx*dx + dy*dy))
    return coords, dist


# ── Formulation 1: DFJ with lazy callbacks (integer solutions only) ────────────
def solve_dfj_lazy(n, dist):
    """Undirected TSP: subtour elimination added lazily on integer solutions."""
    m = gp.Model()
    m.Params.OutputFlag = 0
    m.Params.TimeLimit = 120

    x = {}
    for i in range(n):
        for j in range(i+1, n):
            x[i,j] = m.addVar(vtype=GRB.BINARY, obj=dist[i][j], name=f"x({i},{j})")

    for i in range(n):
        edges = [x[min(i,j), max(i,j)] for j in range(n) if j != i]
        m.addConstr(gp.quicksum(edges) == 2, name=f"deg({i})")

    def subtour_elim_lazy(model, where):
        if where != GRB.Callback.MIPSOL:
            return
        vals = model.cbGetSolution(x)
        adj = [[] for _ in range(n)]
        for (i,j), v in vals.items():
            if v > 0.5:
                adj[i].append(j); adj[j].append(i)
        visited = [False]*n
        for start in range(n):
            if visited[start]: continue
            comp, stack = [], [start]
            while stack:
                node = stack.pop()
                if visited[node]: continue
                visited[node] = True
                comp.append(node)
                stack.extend(adj[node])
            if len(comp) < n:
                edges_in = [x[min(i,j),max(i,j)] for i in comp for j in comp if i < j]
                model.cbLazy(gp.quicksum(edges_in) <= len(comp)-1)

    m.Params.LazyConstraints = 1
    m.optimize(subtour_elim_lazy)
    assert m.Status == GRB.OPTIMAL, f"DFJ-lazy: not optimal (status={m.Status})"
    return int(round(m.ObjVal))


# ── Formulation 2: DFJ with user cuts on LP fractional solutions ───────────────
def _find_subtours(n, xvals):
    """BFS connected components on the fractional solution graph."""
    adj = [[] for _ in range(n)]
    for (i,j), v in xvals.items():
        if v > 1e-6:
            adj[i].append((j, v)); adj[j].append((i, v))
    visited = [False]*n
    comps = []
    for start in range(n):
        if visited[start]: continue
        comp, stack = [], [start]
        while stack:
            node = stack.pop()
            if visited[node]: continue
            visited[node] = True
            comp.append(node)
            stack.extend(nb for nb,_ in adj[node])
        comps.append(comp)
    return comps

def solve_dfj_cuts(n, dist):
    """Undirected TSP: subtour elimination via user cuts on LP relaxation."""
    m = gp.Model()
    m.Params.OutputFlag = 0
    m.Params.TimeLimit = 120
    m.Params.PreCrush = 1  # required when adding user cuts

    x = {}
    for i in range(n):
        for j in range(i+1, n):
            x[i,j] = m.addVar(vtype=GRB.BINARY, obj=dist[i][j], name=f"x({i},{j})")

    for i in range(n):
        edges = [x[min(i,j), max(i,j)] for j in range(n) if j != i]
        m.addConstr(gp.quicksum(edges) == 2, name=f"deg({i})")

    def subtour_elim_cuts(model, where):
        # Add SECs as user cuts on fractional LP solutions
        if where == GRB.Callback.MIPNODE and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
            vals = model.cbGetNodeRel(x)
            comps = _find_subtours(n, vals)
            for comp in comps:
                if len(comp) < n:
                    edges_in = [x[min(i,j),max(i,j)] for i in comp for j in comp if i < j]
                    model.cbCut(gp.quicksum(edges_in) <= len(comp)-1)
        # Also add lazy cuts on integer solutions
        if where == GRB.Callback.MIPSOL:
            vals = model.cbGetSolution(x)
            comps = _find_subtours(n, vals)
            for comp in comps:
                if len(comp) < n:
                    edges_in = [x[min(i,j),max(i,j)] for i in comp for j in comp if i < j]
                    model.cbLazy(gp.quicksum(edges_in) <= len(comp)-1)

    m.Params.LazyConstraints = 1
    m.optimize(subtour_elim_cuts)
    assert m.Status == GRB.OPTIMAL, f"DFJ-cuts: not optimal (status={m.Status})"
    return int(round(m.ObjVal))


# ── Formulation 3: MTZ compact (directed arcs, no callbacks) ──────────────────
def solve_mtz(n, dist):
    """Directed TSP: Miller-Tucker-Zemlin position variables (compact MIP)."""
    m = gp.Model()
    m.Params.OutputFlag = 0
    m.Params.TimeLimit = 120

    # y[i,j] = 1 if directed arc i->j is in tour
    y = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                y[i,j] = m.addVar(vtype=GRB.BINARY, obj=dist[i][j], name=f"y({i},{j})")

    # Position variables u[i], i=1..n-1 (node 0 is depot)
    u = {}
    for i in range(1, n):
        u[i] = m.addVar(lb=1, ub=n-1, vtype=GRB.CONTINUOUS, name=f"u({i})")

    # Out-degree == 1
    for i in range(n):
        m.addConstr(gp.quicksum(y[i,j] for j in range(n) if j != i) == 1, name=f"out({i})")
    # In-degree == 1
    for j in range(n):
        m.addConstr(gp.quicksum(y[i,j] for i in range(n) if i != j) == 1, name=f"in({j})")

    # MTZ subtour-elimination constraints
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                m.addConstr(u[i] - u[j] + n * y[i,j] <= n-1, name=f"mtz({i},{j})")

    m.optimize()
    assert m.Status == GRB.OPTIMAL, f"MTZ: not optimal (status={m.Status})"
    return int(round(m.ObjVal))


# ── Solve all instances with all three formulations ────────────────────────────
INSTANCES = [
    ("n08_s1", 8,  42),
    ("n10_s1", 10, 42),
    ("n12_s1", 12, 42),
    ("n15_s1", 15, 42),
    ("n20_s1", 20, 42),
    ("n08_s2", 8,  137),
    ("n12_s2", 12, 137),
    ("n15_s2", 15, 137),
]

print("Solving TSP instances with Gurobi (3 formulations for cross-validation)...")
results = []
all_ok = True
for name, n, seed in INSTANCES:
    coords, dist = make_instance(n, seed)

    opt_lazy = solve_dfj_lazy(n, dist)
    opt_cuts = solve_dfj_cuts(n, dist)
    opt_mtz  = solve_mtz(n, dist)

    agree = (opt_lazy == opt_cuts == opt_mtz)
    status = "OK" if agree else "MISMATCH"
    print(f"  {name}: n={n} seed={seed}  DFJ-lazy={opt_lazy}  DFJ-cuts={opt_cuts}  MTZ={opt_mtz}  [{status}]")

    if not agree:
        all_ok = False
    else:
        results.append((name, n, seed, opt_lazy, coords, dist))

if not all_ok:
    print("\nERROR: formulation disagreement — fix instances before regenerating fixtures.")
    sys.exit(1)

print(f"\nAll {len(results)} instances cross-validated OK.")

# ---- Generate C source ----
lines = []
lines.append("/* Auto-generated TSP test fixtures. Do not edit manually.")
lines.append("   Generated by test/gen_tsp_fixtures.py using Gurobi 13.")
lines.append("   Instances: random Euclidean TSP, rounded integer distances.")
lines.append("*/")
lines.append("")

lines.append("/* ---- Instance definitions ---- */")
lines.append("typedef struct { int n; int seed; int opt; const int *dist; } TspInstance;")
lines.append("")

for name, n, seed, opt, coords, dist in results:
    lines.append(f"/* {name}: n={n} seed={seed} opt={opt} */")
    flat = []
    for row in dist:
        flat.extend(row)
    arr = ", ".join(str(v) for v in flat)
    lines.append(f"static const int dist_{name}[{n*n}] = {{{arr}}};")
    lines.append("")

lines.append("static const TspInstance TSP_INSTANCES[] = {")
for name, n, seed, opt, coords, dist in results:
    lines.append(f'  {{ {n}, {seed}, {opt}, dist_{name} }},')
lines.append("};")
lines.append(f"#define N_TSP_INSTANCES {len(results)}")

print("\n--- C fixture header ---")
print("\n".join(lines))

# Save
with open("/tmp/tsp_fixtures.h", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\nSaved to /tmp/tsp_fixtures.h")
