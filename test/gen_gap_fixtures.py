#!/usr/bin/env python3
"""
Generate Generalized Assignment Problem (GAP) fixtures.

Problem:
  minimize   sum_{i,j} c[i][j] * x[i][j]
  subject to sum_i x[i][j] = 1              for all j  (each task assigned once)
             sum_j w[i][j] * x[i][j] <= b[i] for all i  (agent capacity)
             x[i][j] in {0,1}

where:
  - c[i][j] = cost of assigning task j to agent i
  - w[i][j] = resource use of task j at agent i
  - b[i]    = capacity of agent i

Diversity axes:
  - agents/tasks ratio: few-agents/many-tasks (hard) vs balanced
  - capacity tightness: tight (alpha~1.2) vs medium (alpha~1.5) vs loose (alpha~1.9)
  - cost type: uniform, structured (agent specialization), correlated (c~w),
               anti-correlated (c~1/w), skewed (one cheap agent)
  - weight type: uniform [1,10], variable [1,20], heavy [5,20]

Naming: gap_{A}a{T}t_sd{seed}_{tightness}_{cost_type}_{weight_type}
  A = number of agents, T = number of tasks
"""

import sys
import math
import random
import gzip
from pathlib import Path


def generate_gap_instance(n_agents, n_tasks, seed,
                          alpha,
                          cost_type="uniform",
                          weight_type="uniform"):
    """
    Generate a GAP instance.

    alpha: capacity tightness factor.
      b[i] = round(alpha * sum_j(w[i][j]) / n_agents)
      With alpha=1.2 the capacity is ~20% above the fair-share average demand
      (tight), with alpha=1.5 it is 50% above (medium).

    cost_type:
      "uniform"       — c[i][j] ~ Uniform[1,100], independent of weights
      "correlated"    — c[i][j] proportional to w[i][j] (heavy=expensive)
      "anticorrelated"— c[i][j] inversely proportional (heavy=cheap)
      "structured"    — agent i has low cost for its preferred task group
      "skewed"        — agent 0 is cheapest for most tasks; others are expensive

    weight_type:
      "uniform"  — w[i][j] ~ Uniform[1,10]
      "variable" — w[i][j] ~ Uniform[1,20]  (wider spread, harder packing)
      "heavy"    — w[i][j] ~ Uniform[5,20]  (all tasks are heavy)
    """
    rng = random.Random(seed)

    # --- weights ---
    if weight_type == "uniform":
        w_lo, w_hi = 1, 10
    elif weight_type == "variable":
        w_lo, w_hi = 1, 20
    elif weight_type == "heavy":
        w_lo, w_hi = 5, 20
    else:
        raise ValueError(f"Unknown weight_type: {weight_type}")

    weights = [[rng.randint(w_lo, w_hi) for _ in range(n_tasks)]
               for _ in range(n_agents)]

    # --- capacities ---
    # b[i] = round(alpha * sum_j w[i][j] / n_agents)
    # Each agent i owns the sum of its row; dividing by n_agents gives the
    # "fair share" if tasks were split evenly across agents using agent i's
    # specific weights.  Multiplying by alpha adds slack (alpha>1 = feasible).
    capacity = [
        max(1, round(alpha * sum(weights[i]) / n_agents))
        for i in range(n_agents)
    ]

    # --- costs ---
    costs = [[0] * n_tasks for _ in range(n_agents)]

    if cost_type == "uniform":
        for i in range(n_agents):
            for j in range(n_tasks):
                costs[i][j] = rng.randint(1, 100)

    elif cost_type == "correlated":
        # c[i][j] = w[i][j] * scale + small noise  →  heavy tasks expensive
        for i in range(n_agents):
            for j in range(n_tasks):
                scale = rng.uniform(3, 7)
                noise = rng.randint(0, 10)
                costs[i][j] = max(1, int(weights[i][j] * scale + noise))

    elif cost_type in ("anticorrelated", "anticorr"):
        # c[i][j] inversely proportional to w[i][j]: light tasks are expensive.
        # Creates the largest LP-MIP gap: LP can fractionally use cheap (heavy)
        # tasks for every agent, but integrality + capacity force harder choices.
        for i in range(n_agents):
            for j in range(n_tasks):
                inv = (w_hi + w_lo + 1 - weights[i][j])
                costs[i][j] = max(1, rng.randint(1, 20) + inv * 5)

    elif cost_type == "structured":
        # Divide tasks into n_agents groups; agent i has low cost for group i.
        # This creates natural specialization — optimal assignment exploits it.
        group_size = math.ceil(n_tasks / n_agents)
        for i in range(n_agents):
            preferred_lo = i * group_size
            preferred_hi = min((i + 1) * group_size, n_tasks)
            for j in range(n_tasks):
                if preferred_lo <= j < preferred_hi:
                    costs[i][j] = rng.randint(1, 20)   # cheap for own group
                else:
                    costs[i][j] = rng.randint(60, 120)  # expensive otherwise

    elif cost_type == "skewed":
        # Agent 0 has uniformly low costs; agents 1..n-1 are expensive.
        # Creates high LP/MIP gap because fractional assignment to agent 0 is
        # optimal in LP but integrality + capacity make it infeasible in MIP.
        for i in range(n_agents):
            for j in range(n_tasks):
                if i == 0:
                    costs[i][j] = rng.randint(1, 15)
                else:
                    costs[i][j] = rng.randint(50, 120)
    else:
        raise ValueError(f"Unknown cost_type: {cost_type}")

    return costs, weights, capacity


def write_gap_mps(filepath, n_agents, n_tasks, costs, weights, capacity):
    """Write a GAP instance to an MPS text file."""
    with open(filepath, 'w') as f:
        f.write("NAME          GAP\n")

        f.write("ROWS\n")
        f.write(" N  OBJ\n")
        for j in range(n_tasks):
            f.write(f" E  ASGN{j}\n")
        for i in range(n_agents):
            f.write(f" L  CAP{i}\n")

        f.write("COLUMNS\n")
        for i in range(n_agents):
            for j in range(n_tasks):
                var = f"x{i}_{j}"
                f.write(f"    {var:<12s}  OBJ    {costs[i][j]}\n")
                f.write(f"    {var:<12s}  ASGN{j:<6d}  1\n")
                f.write(f"    {var:<12s}  CAP{i:<7d}  {weights[i][j]}\n")

        f.write("RHS\n")
        for j in range(n_tasks):
            f.write(f"    RHS       ASGN{j}  1\n")
        for i in range(n_agents):
            f.write(f"    RHS       CAP{i}  {capacity[i]}\n")

        f.write("BOUNDS\n")
        for i in range(n_agents):
            for j in range(n_tasks):
                f.write(f" BV BND       x{i}_{j}\n")

        f.write("ENDATA\n")


def solve_gap_gurobi(n_agents, n_tasks, costs, weights, capacity, timeout=60):
    """Solve GAP with Gurobi to obtain certified optimum."""
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError:
        print("  WARNING: gurobipy not available", file=sys.stderr)
        return None, "no_solver"

    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.setParam("TimeLimit", timeout)
    env.start()

    m = gp.Model(env=env)

    x = [[m.addVar(vtype=GRB.BINARY, name=f"x{i}_{j}")
          for j in range(n_tasks)]
         for i in range(n_agents)]

    m.setObjective(
        gp.quicksum(costs[i][j] * x[i][j]
                    for i in range(n_agents)
                    for j in range(n_tasks)),
        GRB.MINIMIZE
    )

    # Assignment: each task to exactly one agent
    for j in range(n_tasks):
        m.addConstr(gp.quicksum(x[i][j] for i in range(n_agents)) == 1,
                    name=f"assign_{j}")

    # Capacity per agent
    for i in range(n_agents):
        m.addConstr(
            gp.quicksum(weights[i][j] * x[i][j] for j in range(n_tasks))
            <= capacity[i],
            name=f"cap_{i}"
        )

    m.optimize()

    status_map = {
        GRB.OPTIMAL:    "optimal",
        GRB.INFEASIBLE: "infeasible",
        GRB.TIME_LIMIT: "timelimit",
    }
    status = status_map.get(m.Status, f"gurobi_{m.Status}")

    obj = None
    if m.SolCount > 0:
        obj = round(m.ObjVal)

    env.dispose()
    return obj, status


def tightness_label(alpha):
    if alpha < 1.3:
        return "tight"
    elif alpha < 1.65:
        return "medium"
    else:
        return "loose"


# ---------------------------------------------------------------------------
# Instance specifications
# (n_agents, n_tasks, seed, alpha, cost_type, weight_type, gurobi_timeout)
#
# Diversity axes covered:
#   cost_type:
#     anticorr   — cheap tasks are heavy; largest LP-MIP gap; requires branching
#     uniform    — random independent costs; tests basic correctness
#     structured — agent specialisation; tests assignment-side cuts
#     skewed     — one cheap agent; tests preprocessing / fixing
#     correlated — heavy tasks expensive; cost-weight interaction
#   weight_type:
#     uniform [1,10], variable [1,20], heavy [5,20]
#   capacity tightness:
#     tight  alpha≈1.05-1.15  medium alpha≈1.3-1.4  loose alpha≈1.6+
# ---------------------------------------------------------------------------
SPECS = [
    # ── Anti-correlated: require genuine B&B ─────────────────────────────────
    # Warm-up: smallest, needs ~500 nodes
    (4, 35, 42,  1.08, "anticorr", "uniform",  30),
    # Small-medium: ~1000-5000 nodes, < 1s
    (6, 30, 42,  1.15, "anticorr", "uniform",  30),
    (4, 35, 137, 1.08, "anticorr", "uniform",  30),
    # Medium: ~5000-12000 nodes, 1-3s
    (4, 40, 42,  1.10, "anticorr", "variable", 60),
    (6, 35, 42,  1.15, "anticorr", "uniform",  60),
    (5, 30, 42,  1.12, "anticorr", "variable", 60),
    # Hard: ~10000-25000 nodes, 3-8s; heavy weights make packing harder
    (6, 30, 137, 1.15, "anticorr", "heavy",    90),
    (5, 35, 137, 1.12, "anticorr", "variable", 90),

    # ── Other cost types: verify correctness (solve at root) ─────────────────
    # Tests that uniform-cost GAP is handled correctly (no spurious infeasibility)
    (5, 20, 42,  1.20, "uniform",    "uniform",  30),
    # Tests agent-specialisation structure (structured costs)
    (5, 20, 137, 1.20, "structured", "uniform",  30),
    # Tests the larger, balanced case (many agents, medium slack)
    (10, 30, 42, 1.50, "uniform",    "uniform",  30),
]


def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    results = []

    for (n_ag, n_tk, seed, alpha, cost_type, weight_type, timeout) in SPECS:
        tag = tightness_label(alpha)
        name = f"gap_{n_ag}a{n_tk}t_sd{seed}_{tag}_{cost_type}_{weight_type}"
        print(f"Generating {name} …", flush=True)

        costs, weights, capacity = generate_gap_instance(
            n_ag, n_tk, seed, alpha, cost_type, weight_type
        )

        obj, status = solve_gap_gurobi(n_ag, n_tk, costs, weights, capacity, timeout)

        if status == "optimal" and obj is not None:
            n_vars = n_ag * n_tk
            print(f"  obj={obj}  status=optimal  vars={n_vars} ({n_ag}×{n_tk})")

            mps_path = f"/tmp/{name}.mps"
            write_gap_mps(mps_path, n_ag, n_tk, costs, weights, capacity)

            gz_path = fixture_dir / f"{name}.mps.gz"
            with open(mps_path, 'rb') as fin, gzip.open(gz_path, 'wb') as fout:
                fout.writelines(fin)

            print(f"  → {gz_path}")
            results.append((name, n_ag, n_tk, obj, timeout))

        elif status == "infeasible":
            print(f"  INFEASIBLE — skipping (alpha={alpha} too tight?)")
        else:
            print(f"  {status} obj={obj} — skipping")

    # ── Print C snippet ──────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("Paste into test/CInterfaceTest_gap.c  (GapTestCase array):")
    print("=" * 72)
    print("static const GapTestCase gap_test_cases[] = {")
    for (name, n_ag, n_tk, obj, timeout) in results:
        print(f'  {{"{name}", {obj}, 50000, {timeout}}},')
    print("};")


if __name__ == "__main__":
    main()
