#!/usr/bin/env python3
"""Generate Fixed Charge Network Flow (FCNF) test fixtures with certified optima.

FCNF formulation:
  min  sum_{ij} (c_ij * x_ij + f_ij * y_ij)
  s.t. sum_j x_ji - sum_j x_ij = b_i           (flow balance)
       0 <= x_ij <= u_ij * y_ij                (capacity + linking)
       y_ij in {0,1}, x_ij continuous
  where b_i > 0 = supply, b_i < 0 = demand, b_i = 0 = transshipment.

The x <= u*y linking constraint is the variable-upper-bound (VUB) structure
that exercises flow cuts, knapsack-with-VUB, and MIR generators.

Diversity axes:
  - Topology: grid, layered (multi-stage), hub-and-spoke, complete bipartite,
    Erdos-Renyi random.
  - Demand pattern: single-source/single-sink, distribution, transportation.
  - Cost regime: flat (c~f), fixed-dominated (f>>c), variable-dominated (c>>f),
    heterogeneous magnitudes within one instance.
  - Capacity scale: integer (u in {1..10}), continuous (u in [1, 50]),
    mixed (some arcs u~5, others u~500).

Each fixture is solved by an external solver to certify the optimum.
The certified value is recorded and used by the test only when MIPster
proves optimality (otherwise feasibility alone is required).

Sizing keeps each instance solvable on the dev box in seconds; the medium
random instance is intentionally hard enough to exhaust its node budget,
exercising the feasible-not-proven pass path.
"""

import gzip
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path


# ----------------------------------------------------------------------
# Topology builders. Each returns (n_nodes, arcs, supplies) where:
#   arcs     = list of (tail, head, c, f, u)
#   supplies = dict node -> b_i (supply > 0, demand < 0)
# ----------------------------------------------------------------------


def _scale_supplies_to_capacity(arcs, supplies, slack=0.5):
    """Rescale supplies/demands so total supply <= slack * min(out-cap of sources,
    in-cap of sinks). Ensures the instance is feasible with room for routing.

    `slack` is the fraction of the limiting capacity to use as total flow:
      slack=0.5 means total supply = 50% of the limit (loose).
      slack=0.85 means tight.
    Preserves the relative ratios between sources and between sinks."""
    out_cap = {}
    in_cap = {}
    for (i, j, _c, _f, u) in arcs:
        out_cap[i] = out_cap.get(i, 0.0) + u
        in_cap[j] = in_cap.get(j, 0.0) + u

    sources = [(n, b) for n, b in supplies.items() if b > 0]
    sinks = [(n, -b) for n, b in supplies.items() if b < 0]
    if not sources or not sinks:
        return supplies

    # Total flow F = scale * total_supply, must satisfy
    #   for each source n: scale * b_n <= out_cap[n]   =>  scale <= out_cap[n]/b_n
    #   for each sink n:   scale * b_n <= in_cap[n]    =>  scale <= in_cap[n]/b_n
    scale_max = math.inf
    for n, b in sources:
        if out_cap.get(n, 0.0) > 0:
            scale_max = min(scale_max, out_cap[n] / b)
        else:
            return None  # source has no outgoing arcs — instance broken
    for n, b in sinks:
        if in_cap.get(n, 0.0) > 0:
            scale_max = min(scale_max, in_cap[n] / b)
        else:
            return None  # sink has no incoming arcs

    scale = scale_max * slack
    new_supplies = {}
    for n, b in sources:
        new_supplies[n] = max(1, int(round(b * scale)))
    for n, b in sinks:
        new_supplies[n] = -max(1, int(round(b * scale)))

    # Re-balance so total supply == total demand. Adjust the largest source.
    total_in = sum(v for v in new_supplies.values() if v > 0)
    total_out = -sum(v for v in new_supplies.values() if v < 0)
    diff = total_in - total_out
    if diff != 0:
        # Adjust the source with most slack (largest supply)
        src_n = max(new_supplies, key=lambda k: new_supplies[k])
        new_supplies[src_n] -= diff
        if new_supplies[src_n] <= 0:
            # Try absorbing into a sink instead
            new_supplies[src_n] += diff  # undo
            sink_n = min(new_supplies, key=lambda k: new_supplies[k])
            new_supplies[sink_n] += diff
            if new_supplies[sink_n] >= 0:
                return None
    # Sanity: re-verify capacity constraints with new integer values
    for n, b in [(n, b) for n, b in new_supplies.items() if b > 0]:
        if b > out_cap.get(n, 0.0):
            return None
    for n, b in [(n, -b) for n, b in new_supplies.items() if b < 0]:
        if b > in_cap.get(n, 0.0):
            return None
    return new_supplies

def build_grid(rows, cols, supply_nodes, demand_nodes, rng,
               c_range, f_range, u_range, u_continuous=False):
    """Rectangular grid; arcs only between 4-neighbors, both directions."""
    n_nodes = rows * cols
    arcs = []
    seen = set()
    for r in range(rows):
        for c in range(cols):
            i = r * cols + c
            for (dr, dc) in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                rr, cc = r + dr, c + dc
                if 0 <= rr < rows and 0 <= cc < cols:
                    j = rr * cols + cc
                    if (i, j) in seen:
                        continue
                    seen.add((i, j))
                    cc_cost = rng.uniform(*c_range)
                    ff = rng.uniform(*f_range)
                    if u_continuous:
                        uu = rng.uniform(*u_range)
                    else:
                        uu = float(rng.randint(int(u_range[0]), int(u_range[1])))
                    arcs.append((i, j, round(cc_cost, 2), round(ff, 2),
                                 round(uu, 2)))
    supplies = {}
    for n, b in supply_nodes:
        supplies[n] = supplies.get(n, 0) + b
    for n, b in demand_nodes:
        supplies[n] = supplies.get(n, 0) - b
    return n_nodes, arcs, supplies


def build_layered(layer_sizes, rng, c_range, f_range, u_range,
                  u_continuous=False, sources_per_first=2, sinks_per_last=2):
    """Multi-stage transshipment: nodes grouped in layers, arcs only forward."""
    n_nodes = sum(layer_sizes)
    starts = [0]
    for s in layer_sizes[:-1]:
        starts.append(starts[-1] + s)

    arcs = []
    for li in range(len(layer_sizes) - 1):
        sz1, sz2 = layer_sizes[li], layer_sizes[li + 1]
        for a in range(sz1):
            i = starts[li] + a
            for b in range(sz2):
                j = starts[li + 1] + b
                cc = rng.uniform(*c_range)
                ff = rng.uniform(*f_range)
                if u_continuous:
                    uu = rng.uniform(*u_range)
                else:
                    uu = float(rng.randint(int(u_range[0]), int(u_range[1])))
                arcs.append((i, j, round(cc, 2), round(ff, 2), round(uu, 2)))

    # First layer = sources, last layer = sinks
    supplies = {}
    first_size = layer_sizes[0]
    last_size = layer_sizes[-1]
    last_start = starts[-1]
    src_indices = list(range(first_size))
    snk_indices = list(range(last_start, last_start + last_size))
    rng.shuffle(src_indices)
    rng.shuffle(snk_indices)

    src_picked = src_indices[:sources_per_first]
    snk_picked = snk_indices[:sinks_per_last]

    # Per-source supply, balanced
    total = 10 * len(snk_picked)
    per_src = total // len(src_picked)
    rem = total - per_src * len(src_picked)
    for k, s in enumerate(src_picked):
        supplies[s] = per_src + (1 if k < rem else 0)
    per_snk = total // len(snk_picked)
    rem2 = total - per_snk * len(snk_picked)
    for k, t in enumerate(snk_picked):
        supplies[t] = -(per_snk + (1 if k < rem2 else 0))
    return n_nodes, arcs, supplies


def build_hub_spoke(n_total, n_hubs, rng, c_range, f_range, u_range,
                    u_continuous=False, hub_arc_capacity_mult=10):
    """Scale-free-ish: a few high-degree hubs connected to all spokes.
    Every spoke connects to every hub (both directions).
    Hubs interconnect (full mesh, both directions).
    Arcs through hubs have larger capacity to make hub structure useful."""
    arcs = []
    hubs = list(range(n_hubs))
    spokes = list(range(n_hubs, n_total))

    # Spoke -> hub and hub -> spoke
    for h in hubs:
        for s in spokes:
            for (i, j) in [(s, h), (h, s)]:
                cc = rng.uniform(*c_range)
                ff = rng.uniform(*f_range)
                if u_continuous:
                    uu = rng.uniform(u_range[0] * hub_arc_capacity_mult,
                                     u_range[1] * hub_arc_capacity_mult)
                else:
                    uu = float(rng.randint(
                        int(u_range[0] * hub_arc_capacity_mult),
                        int(u_range[1] * hub_arc_capacity_mult)))
                arcs.append((i, j, round(cc, 2), round(ff, 2), round(uu, 2)))

    # Hub <-> hub
    for a in range(n_hubs):
        for b in range(n_hubs):
            if a == b:
                continue
            cc = rng.uniform(*c_range)
            ff = rng.uniform(*f_range)
            if u_continuous:
                uu = rng.uniform(u_range[0] * hub_arc_capacity_mult,
                                 u_range[1] * hub_arc_capacity_mult)
            else:
                uu = float(rng.randint(
                    int(u_range[0] * hub_arc_capacity_mult),
                    int(u_range[1] * hub_arc_capacity_mult)))
            arcs.append((hubs[a], hubs[b],
                         round(cc, 2), round(ff, 2), round(uu, 2)))

    # Demand: all spokes send to one chosen sink spoke
    sink = spokes[-1]
    sources = spokes[:-1]
    supplies = {sink: -len(sources) * 3}
    for s in sources:
        supplies[s] = 3
    return n_total, arcs, supplies


def build_complete_bipartite(n_left, n_right, rng,
                             c_range, f_range, u_range, u_continuous=False):
    """Transportation-style: left = sources, right = sinks, full bipartite arcs."""
    n_nodes = n_left + n_right
    arcs = []
    for i in range(n_left):
        for j in range(n_right):
            cc = rng.uniform(*c_range)
            ff = rng.uniform(*f_range)
            if u_continuous:
                uu = rng.uniform(*u_range)
            else:
                uu = float(rng.randint(int(u_range[0]), int(u_range[1])))
            arcs.append((i, n_left + j, round(cc, 2), round(ff, 2),
                         round(uu, 2)))

    # Total supply == total demand
    supplies_left = []
    total = 5 * n_right
    per = total // n_left
    rem = total - per * n_left
    for k in range(n_left):
        supplies_left.append(per + (1 if k < rem else 0))
    supplies_right = []
    per_r = total // n_right
    rem_r = total - per_r * n_right
    for k in range(n_right):
        supplies_right.append(per_r + (1 if k < rem_r else 0))

    supplies = {}
    for i in range(n_left):
        supplies[i] = supplies_left[i]
    for j in range(n_right):
        supplies[n_left + j] = -supplies_right[j]
    return n_nodes, arcs, supplies


def build_random(n_nodes, avg_out_degree, rng,
                 c_range, f_range, u_range, u_continuous=False,
                 n_sources=2, n_sinks=3, capacity_tightness=1.0):
    """Erdos-Renyi-style random digraph. capacity_tightness in (0, 1]:
    1.0 = loose, 0.5 = each arc capacity halved => total cap barely covers demand."""
    arcs = []
    seen = set()
    target_arcs = n_nodes * avg_out_degree
    attempts = 0
    while len(arcs) < target_arcs and attempts < target_arcs * 20:
        i = rng.randrange(n_nodes)
        j = rng.randrange(n_nodes)
        if i == j or (i, j) in seen:
            attempts += 1
            continue
        seen.add((i, j))
        cc = rng.uniform(*c_range)
        ff = rng.uniform(*f_range)
        if u_continuous:
            uu = rng.uniform(*u_range) * capacity_tightness
        else:
            uu_int = rng.randint(int(u_range[0]), int(u_range[1]))
            uu = max(1.0, uu_int * capacity_tightness)
        arcs.append((i, j, round(cc, 2), round(ff, 2), round(uu, 2)))
        attempts += 1

    indices = list(range(n_nodes))
    rng.shuffle(indices)
    src = indices[:n_sources]
    snk = indices[n_sources:n_sources + n_sinks]

    supplies = {}
    total = 8 * n_sinks
    per_s = total // n_sources
    rem = total - per_s * n_sources
    for k, s in enumerate(src):
        supplies[s] = per_s + (1 if k < rem else 0)
    per_t = total // n_sinks
    rem2 = total - per_t * n_sinks
    for k, t in enumerate(snk):
        supplies[t] = -(per_t + (1 if k < rem2 else 0))
    return n_nodes, arcs, supplies


def build_random_heterogeneous(n_nodes, avg_out_degree, rng,
                               n_sources=2, n_sinks=3):
    """Random topology with explicitly heterogeneous magnitudes:
    some arcs have c~O(1), f~O(100); others have c~O(50), f~O(2)."""
    arcs = []
    seen = set()
    target_arcs = n_nodes * avg_out_degree
    attempts = 0
    while len(arcs) < target_arcs and attempts < target_arcs * 20:
        i = rng.randrange(n_nodes)
        j = rng.randrange(n_nodes)
        if i == j or (i, j) in seen:
            attempts += 1
            continue
        seen.add((i, j))
        # Heterogeneous: half "fixed-heavy", half "variable-heavy"
        if rng.random() < 0.5:
            cc = rng.uniform(0.5, 3.0)
            ff = rng.uniform(50.0, 200.0)
            uu = float(rng.randint(2, 8))
        else:
            cc = rng.uniform(20.0, 80.0)
            ff = rng.uniform(1.0, 5.0)
            uu = float(rng.randint(2, 8))
        arcs.append((i, j, round(cc, 2), round(ff, 2), round(uu, 2)))
        attempts += 1

    indices = list(range(n_nodes))
    rng.shuffle(indices)
    src = indices[:n_sources]
    snk = indices[n_sources:n_sources + n_sinks]

    supplies = {}
    total = 6 * n_sinks
    per_s = total // n_sources
    rem = total - per_s * n_sources
    for k, s in enumerate(src):
        supplies[s] = per_s + (1 if k < rem else 0)
    per_t = total // n_sinks
    rem2 = total - per_t * n_sinks
    for k, t in enumerate(snk):
        supplies[t] = -(per_t + (1 if k < rem2 else 0))
    return n_nodes, arcs, supplies


def build_grid_multimagnitude(rows, cols, rng):
    """4x4 grid where some arcs have u~5 (small) and others u~500 (large).
    Tests numerical robustness of x <= u*y inside the LP relaxation."""
    n_nodes = rows * cols
    arcs = []
    seen = set()
    for r in range(rows):
        for c in range(cols):
            i = r * cols + c
            for (dr, dc) in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                rr, cc_ = r + dr, c + dc
                if 0 <= rr < rows and 0 <= cc_ < cols:
                    j = rr * cols + cc_
                    if (i, j) in seen:
                        continue
                    seen.add((i, j))
                    if rng.random() < 0.3:
                        uu = float(rng.randint(200, 800))
                        cc_cost = rng.uniform(0.05, 0.2)
                        ff = rng.uniform(20.0, 50.0)
                    else:
                        uu = float(rng.randint(3, 8))
                        cc_cost = rng.uniform(1.0, 5.0)
                        ff = rng.uniform(5.0, 20.0)
                    arcs.append((i, j, round(cc_cost, 2),
                                 round(ff, 2), round(uu, 2)))

    src = 0
    snk = n_nodes - 1
    supplies = {src: 15, snk: -15}
    return n_nodes, arcs, supplies


# ----------------------------------------------------------------------
# MPS writer
# ----------------------------------------------------------------------

def write_mps(name, n_nodes, arcs, supplies, output_path):
    """Write FCNF as MPS. Variables: x_a (continuous), y_a (binary) for each arc.
    Constraints:
      flow_i: sum(in) - sum(out) = b_i      for each node
      link_a: x_a - u_a * y_a <= 0          for each arc
    """
    n_arcs = len(arcs)

    var_x = [f"x_{a}" for a in range(n_arcs)]
    var_y = [f"y_{a}" for a in range(n_arcs)]
    row_flow = [f"flow_{i}" for i in range(n_nodes)]
    row_link = [f"link_{a}" for a in range(n_arcs)]

    # Build per-variable column entries
    x_entries = [[] for _ in range(n_arcs)]
    y_entries = [[] for _ in range(n_arcs)]
    obj_x = [0.0] * n_arcs
    obj_y = [0.0] * n_arcs
    for a, (tail, head, c, f, u) in enumerate(arcs):
        obj_x[a] = c
        obj_y[a] = f
        # flow balance with convention "outflow - inflow = b_i"
        # (so b_i > 0 = supply, b_i < 0 = demand). At tail x adds +1; at head x adds -1.
        x_entries[a].append((row_flow[tail], 1.0))
        x_entries[a].append((row_flow[head], -1.0))
        # link: x - u*y <= 0
        x_entries[a].append((row_link[a], 1.0))
        y_entries[a].append((row_link[a], -float(u)))

    rhs_flow = [supplies.get(i, 0) for i in range(n_nodes)]

    with gzip.open(output_path, "wt") as fp:
        fp.write(f"NAME          {name}\n")

        fp.write("ROWS\n")
        fp.write(" N  OBJ\n")
        for r in row_flow:
            fp.write(f" E  {r}\n")
        for r in row_link:
            fp.write(f" L  {r}\n")

        fp.write("COLUMNS\n")
        # Continuous x first
        for a in range(n_arcs):
            if obj_x[a] != 0.0:
                fp.write(f"    {var_x[a]:10s}  OBJ       {obj_x[a]}\n")
            for (rname, coef) in x_entries[a]:
                fp.write(f"    {var_x[a]:10s}  {rname:10s}  {coef}\n")
        # Binary y inside INTORG/INTEND markers
        fp.write("    MARKBIN   'MARKER'                 'INTORG'\n")
        for a in range(n_arcs):
            if obj_y[a] != 0.0:
                fp.write(f"    {var_y[a]:10s}  OBJ       {obj_y[a]}\n")
            for (rname, coef) in y_entries[a]:
                fp.write(f"    {var_y[a]:10s}  {rname:10s}  {coef}\n")
        fp.write("    MARKBIN   'MARKER'                 'INTEND'\n")

        fp.write("RHS\n")
        for i, b in enumerate(rhs_flow):
            if b != 0:
                fp.write(f"    RHS       {row_flow[i]:10s}  {b}\n")

        fp.write("BOUNDS\n")
        for a in range(n_arcs):
            fp.write(f" BV BOUND     {var_y[a]}\n")
            # x default lower 0 (fine), no upper bound (link constraint enforces it)

        fp.write("ENDATA\n")


# ----------------------------------------------------------------------
# External certification
# ----------------------------------------------------------------------

def certify_optimum(mps_path, time_limit=600):
    """Solve to optimality with an external solver, return (obj, status)."""
    cert = shutil.which("gurobi_cl")
    if cert is None:
        # Fall back to system cbc
        cert = shutil.which("cbc")
        if cert is None:
            return None, "no external solver found"
        try:
            res = subprocess.run(
                [cert, str(mps_path), "sec", str(time_limit), "solve"],
                capture_output=True, text=True, timeout=time_limit + 60,
            )
            for line in res.stdout.split("\n"):
                if "Optimal - objective value" in line:
                    return float(line.split("value")[-1].strip()), "optimal"
            return None, "cbc could not certify"
        except Exception as e:
            return None, f"cbc error: {e}"

    try:
        res = subprocess.run(
            [cert, f"TimeLimit={time_limit}", "MIPGap=0", str(mps_path)],
            capture_output=True, text=True, timeout=time_limit + 60,
        )
        obj = None
        is_optimal = False
        is_infeasible = False
        for line in res.stdout.split("\n"):
            line = line.strip()
            if line.startswith("Optimal solution found"):
                is_optimal = True
            if line == "Model is infeasible" or line.startswith("Model is infeasible"):
                is_infeasible = True
            if line.startswith("Best objective"):
                # "Best objective 1.234e+02, best bound ..., gap 0.0000%"
                try:
                    obj_str = line.split("Best objective")[1].split(",")[0].strip()
                    obj = float(obj_str)
                except Exception:
                    pass
        if is_infeasible:
            return None, "INFEASIBLE"
        if is_optimal and obj is not None:
            return obj, "optimal"
        return obj, "not proven optimal"
    except Exception as e:
        return None, f"error: {e}"


# ----------------------------------------------------------------------
# Instance roster
# ----------------------------------------------------------------------

def _add(instances, name, n, arcs, supplies, hint, slack):
    """Rescale supplies to fit arc capacities, then append."""
    scaled = _scale_supplies_to_capacity(arcs, supplies, slack=slack)
    if scaled is None:
        raise RuntimeError(f"{name}: rescaling failed (disconnected source/sink?)")
    instances.append((name, n, arcs, scaled, hint))


def build_all(rng_master_seed=42):
    """Returns list of (name, n_nodes, arcs, supplies, hint_max_nodes)."""
    instances = []

    # 1) Tiny: 4x4 grid, 1 source 3 sinks, flat costs, integer u
    rng = random.Random(rng_master_seed + 1)
    n, arcs, sup = build_grid(
        rows=4, cols=4,
        supply_nodes=[(0, 30)],
        demand_nodes=[(5, 10), (10, 10), (15, 10)],
        rng=rng,
        c_range=(1.0, 5.0), f_range=(2.0, 6.0),
        u_range=(3, 12), u_continuous=False,
    )
    _add(instances, "fcnf_grid_4x4_balanced", n, arcs, sup, 200000, slack=0.5)

    # 2) Small: 5x5 grid, 1 src 1 sink, fixed-heavy, integer u
    rng = random.Random(rng_master_seed + 2)
    n, arcs, sup = build_grid(
        rows=5, cols=5,
        supply_nodes=[(0, 12)],
        demand_nodes=[(24, 12)],
        rng=rng,
        c_range=(0.5, 2.0), f_range=(20.0, 60.0),  # f >> c
        u_range=(3, 8), u_continuous=False,
    )
    _add(instances, "fcnf_grid_5x5_fixedheavy", n, arcs, sup, 50000, slack=0.4)

    # 3) Small: 3-stage layered, flat costs, continuous u
    rng = random.Random(rng_master_seed + 3)
    n, arcs, sup = build_layered(
        layer_sizes=[4, 6, 6, 4], rng=rng,
        c_range=(1.0, 5.0), f_range=(3.0, 8.0),
        u_range=(2.0, 12.0), u_continuous=True,
        sources_per_first=2, sinks_per_last=2,
    )
    _add(instances, "fcnf_layered_3stage", n, arcs, sup, 50000, slack=0.25)

    # 4) Tiny: hub-and-spoke 12 nodes, mixed costs and large hub capacity
    rng = random.Random(rng_master_seed + 4)
    n, arcs, sup = build_hub_spoke(
        n_total=12, n_hubs=2, rng=rng,
        c_range=(1.0, 8.0), f_range=(5.0, 30.0),
        u_range=(2, 6), u_continuous=False,
        hub_arc_capacity_mult=10,
    )
    _add(instances, "fcnf_hub_n12", n, arcs, sup, 200000, slack=0.5)

    # 5) Tiny: complete bipartite 5x5, fixed-heavy
    rng = random.Random(rng_master_seed + 5)
    n, arcs, sup = build_complete_bipartite(
        n_left=5, n_right=5, rng=rng,
        c_range=(0.5, 3.0), f_range=(15.0, 50.0),  # f >> c
        u_range=(3, 10), u_continuous=False,
    )
    _add(instances, "fcnf_complete_bip_5x5", n, arcs, sup, 200000, slack=0.3)

    # 6) Small: random ER n=15, flat, continuous
    rng = random.Random(rng_master_seed + 6)
    n, arcs, sup = build_random(
        n_nodes=15, avg_out_degree=3, rng=rng,
        c_range=(1.0, 6.0), f_range=(3.0, 10.0),
        u_range=(3.0, 15.0), u_continuous=True,
        n_sources=2, n_sinks=3, capacity_tightness=1.0,
    )
    _add(instances, "fcnf_random_n15_d3", n, arcs, sup, 50000, slack=0.3)

    # 7) Medium: random ER n=20, heterogeneous magnitudes, tight capacity
    rng = random.Random(rng_master_seed + 7)
    n, arcs, sup = build_random_heterogeneous(
        n_nodes=20, avg_out_degree=4, rng=rng,
        n_sources=2, n_sinks=3,
    )
    _add(instances, "fcnf_random_n20_d4_tight", n, arcs, sup, 5000, slack=0.85)

    # 8) Small: layered, variable-heavy (c >> f), large continuous capacities
    rng = random.Random(rng_master_seed + 8)
    n, arcs, sup = build_layered(
        layer_sizes=[3, 5, 5, 3], rng=rng,
        c_range=(20.0, 80.0), f_range=(0.5, 3.0),  # c >> f
        u_range=(5.0, 25.0), u_continuous=True,
        sources_per_first=2, sinks_per_last=2,
    )
    _add(instances, "fcnf_layered_continuous_loose", n, arcs, sup, 50000, slack=0.25)

    # 9) Tiny: 4x4 grid with mixed-magnitude capacities (some u~5, others u~500)
    rng = random.Random(rng_master_seed + 9)
    n, arcs, sup = build_grid_multimagnitude(rows=4, cols=4, rng=rng)
    _add(instances, "fcnf_grid_4x4_multimagnitude", n, arcs, sup, 200000, slack=0.5)

    return instances


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    print("=== Generating FCNF Fixtures ===\n")

    rows = []
    instances = build_all()
    for name, n_nodes, arcs, supplies, hint_max_nodes in instances:
        path = fixture_dir / f"{name}.mps.gz"
        write_mps(name, n_nodes, arcs, supplies, path)
        n_arcs = len(arcs)
        print(f"  {name}: n={n_nodes}, arcs={n_arcs}  -> {path}")
        rows.append((name, n_nodes, n_arcs, hint_max_nodes, path))

    print()
    print("=== Certifying optima ===\n")
    certified = []
    for (name, n_nodes, n_arcs, hint, path) in rows:
        obj, status = certify_optimum(path, time_limit=600)
        if obj is not None and status == "optimal":
            print(f"  {name}: certified obj = {obj:.6f}")
            certified.append((name, n_nodes, n_arcs, hint, obj))
        else:
            print(f"  {name}: NOT CERTIFIED ({status}, best={obj})")
            certified.append((name, n_nodes, n_arcs, hint, None))

    print()
    print("=== Test case array (paste into CInterfaceTest_fcnf_fixtures.c) ===\n")
    print("static const FCNFTestCase fcnf_test_cases[] = {")
    for (name, n_nodes, n_arcs, hint, obj) in certified:
        if obj is None:
            print(f'  {{ "{name}", NAN, {hint} }},  /* n={n_nodes} arcs={n_arcs} */')
        else:
            print(f'  {{ "{name}", {obj:.6f}, {hint} }},  /* n={n_nodes} arcs={n_arcs} */')
    print("};")
    return 0


if __name__ == "__main__":
    sys.exit(main())
