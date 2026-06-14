#!/usr/bin/env python3
"""Generate UPMSP-ST (Unrelated Parallel Machine Scheduling with
Sequence-Dependent Setup Times) MPS fixtures for MIPster testing.

Formulation (big-M):
  Variables:
    Cmax         continuous, obj=1 (makespan variant only)
    C[j]         continuous, completion time of job j
    x[j][k]      binary, 1 if job j assigned to machine k
    d[i][j][k]   binary, i < j; 1 if job i precedes job j on machine k

  Constraints (i < j, machine k):
    asgn_j:       sum_k x[j][k] = 1
    lb_j_k:       C[j] - M*x[j][k]                             >= p[j][k] - M
    prA_i_j_k:    C[j]-C[i] - M*d[i,j,k] - M*x[i,k] - M*x[j,k] >= p[j,k]+s[i,j,k] - 3M
    prB_i_j_k:    C[i]-C[j] + M*d[i,j,k] - M*x[i,k] - M*x[j,k] >= p[i,k]+s[j,i,k] - 2M
    cmax_j:       Cmax - C[j]                                   >= 0  (makespan only)

  Objective: minimize Cmax  OR  minimize sum_j w[j]*C[j]

Six diverse instances:
  1. n=8, m=2, integer proc/setup, Cmax          (moderate, asymmetric machines)
  2. n=8, m=3, integer proc/setup, WCT           (diverse weights, 3 machines)
  3. n=9, m=3, integer proc/setup, Cmax          (larger, harder)
  4. n=8, m=3, fractional proc/setup (0.5 grid), Cmax  (non-integer optimal)
  5. n=9, m=2, large setup times, Cmax           (setups dominate sequencing)
  6. n=10, m=3, integer proc/setup, WCT          (largest, hardest WCT)
"""

import gzip, random, subprocess, sys
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────────────
# MPS writer
# ──────────────────────────────────────────────────────────────────────────────

def make_upms_mps(name, n, m, proc, setup, weights=None):
    """Return MPS text for a UPMSP-ST instance.

    proc[j][k]       processing time of job j on machine k
    setup[i][j][k]   setup time from job i to job j on machine k (i != j)
    weights[j]       objective weight of job j  (None → Cmax objective)
    """
    cmax_obj = (weights is None)

    # BigM = n * (max processing time + max setup time)
    max_p = max(proc[j][k] for j in range(n) for k in range(m))
    max_s = max(setup[i][j][k]
                for i in range(n) for j in range(n) for k in range(m)
                if i != j)
    BIG = n * (max_p + max_s)

    # ── column index helpers ──────────────────────────────────────────────────
    # Continuous block: Cmax (opt), C[0..n-1]
    col_cmax = 0 if cmax_obj else None
    col_C    = lambda j: (1 + j) if cmax_obj else j
    n_cont   = (1 + n) if cmax_obj else n

    # Binary block: x[j][k], then d[i][j][k] for i<j
    col_x    = lambda j, k: n_cont + j * m + k
    pairs    = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pair_idx = {(i, j): p for p, (i, j) in enumerate(pairs)}
    col_d    = lambda i, j, k: n_cont + n * m + pair_idx[(i, j)] * m + k

    # ── build rows as list of (name, sense, rhs, {col: coef}) ────────────────
    rows = []

    def add(rname, sense, rhs, coefs):
        rows.append((rname, sense, rhs, coefs))

    # Assignment
    for j in range(n):
        add(f"asgn_{j}", "E", 1.0,
            {col_x(j, k): 1.0 for k in range(m)})

    # Lower bound on C[j] given machine assignment
    for j in range(n):
        for k in range(m):
            add(f"lb_{j}_{k}", "G", proc[j][k] - BIG,
                {col_C(j): 1.0, col_x(j, k): -BIG})

    # Precedence A & B for every ordered pair (i<j) and machine k
    for i, j in pairs:
        for k in range(m):
            d = col_d(i, j, k)
            # A: if i before j on k → C[j] ≥ C[i] + p[j][k] + s[i][j][k]
            add(f"prA_{i}_{j}_{k}", "G",
                proc[j][k] + setup[i][j][k] - 3 * BIG,
                {col_C(j): 1.0, col_C(i): -1.0,
                 d: -BIG, col_x(i, k): -BIG, col_x(j, k): -BIG})
            # B: if j before i on k → C[i] ≥ C[j] + p[i][k] + s[j][i][k]
            add(f"prB_{i}_{j}_{k}", "G",
                proc[i][k] + setup[j][i][k] - 2 * BIG,
                {col_C(i): 1.0, col_C(j): -1.0,
                 d: BIG, col_x(i, k): -BIG, col_x(j, k): -BIG})

    # Makespan linking
    if cmax_obj:
        for j in range(n):
            add(f"cmax_{j}", "G", 0.0,
                {col_cmax: 1.0, col_C(j): -1.0})

    # ── build column→rows inverted index ─────────────────────────────────────
    col_rows = {}   # col_idx → [(row_idx, coef)]
    for ri, (_, _, _, coefs) in enumerate(rows):
        for ci, coef in coefs.items():
            col_rows.setdefault(ci, []).append((ri, coef))

    # Objective coefficients per column
    obj = {}
    if cmax_obj:
        obj[col_cmax] = 1.0
    else:
        for j in range(n):
            obj[col_C(j)] = float(weights[j])

    # ── ordered column list ───────────────────────────────────────────────────
    all_cols = []          # (idx, name, "CONT"|"BIN")
    if cmax_obj:
        all_cols.append((col_cmax, "Cmax", "CONT"))
    for j in range(n):
        all_cols.append((col_C(j), f"C_{j}", "CONT"))
    for j in range(n):
        for k in range(m):
            all_cols.append((col_x(j, k), f"x_{j}_{k}", "BIN"))
    for idx, (i, j) in enumerate(pairs):
        for k in range(m):
            all_cols.append((col_d(i, j, k), f"d_{i}_{j}_{k}", "BIN"))

    # ── write MPS ─────────────────────────────────────────────────────────────
    L = []

    L.append(f"NAME          {name}")
    L.append("ROWS")
    L.append(" N  OBJ")
    for rname, sense, _, _ in rows:
        L.append(f" {sense}  {rname}")

    L.append("COLUMNS")
    in_int = False
    mc = 0
    for ci, cname, ctype in all_cols:
        if ctype == "BIN" and not in_int:
            L.append(f"    MARK{mc:04d}  'MARKER'                 'INTORG'")
            mc += 1; in_int = True
        elif ctype == "CONT" and in_int:
            L.append(f"    MARK{mc:04d}  'MARKER'                 'INTEND'")
            mc += 1; in_int = False

        if ci in obj:
            L.append(f"    {cname:<14s}  OBJ           {obj[ci]}")
        for ri, coef in col_rows.get(ci, []):
            rname = rows[ri][0]
            L.append(f"    {cname:<14s}  {rname:<14s}  {coef}")

    if in_int:
        L.append(f"    MARK{mc:04d}  'MARKER'                 'INTEND'")

    L.append("RHS")
    for rname, _, rhs, _ in rows:
        if rhs != 0.0:
            L.append(f"    RHS           {rname:<14s}  {rhs}")

    L.append("BOUNDS")
    for ci, cname, ctype in all_cols:
        if ctype == "BIN":
            L.append(f" BV BOUND         {cname}")

    L.append("ENDATA")
    return "\n".join(L) + "\n"


# ──────────────────────────────────────────────────────────────────────────────
# Instance generator helpers
# ──────────────────────────────────────────────────────────────────────────────

def gen_proc(rng, n, m, lo, hi, asymmetry=True):
    """Generate processing time matrix.  With asymmetry=True, machine 0 is
    generally slower (higher values) to encourage non-trivial assignment."""
    p = [[rng.randint(lo, hi) for _ in range(m)] for _ in range(n)]
    if asymmetry and m >= 2:
        # Make machine 0 systematically slower (multiply by ~1.5)
        for j in range(n):
            p[j][0] = min(hi, int(p[j][0] * 1.5))
    return p


def gen_proc_frac(rng, n, m, lo, hi):
    """Generate processing times on a 0.5 grid in [lo, hi]."""
    return [[rng.choice([x * 0.5 for x in range(int(lo * 2), int(hi * 2) + 1)])
             for _ in range(m)] for _ in range(n)]


def gen_setup(rng, n, m, lo, hi):
    """Generate setup time tensor setup[i][j][k] for i≠j."""
    s = [[[0] * m for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                for k in range(m):
                    s[i][j][k] = rng.randint(lo, hi)
    return s


def gen_setup_frac(rng, n, m, lo, hi):
    """Setup times on 0.5 grid."""
    s = [[[0.0] * m for _ in range(n)] for _ in range(n)]
    choices = [x * 0.5 for x in range(int(lo * 2), int(hi * 2) + 1)]
    for i in range(n):
        for j in range(n):
            if i != j:
                for k in range(m):
                    s[i][j][k] = rng.choice(choices)
    return s


# ──────────────────────────────────────────────────────────────────────────────
# Verification with external solver
# ──────────────────────────────────────────────────────────────────────────────

def verify_with_solver(mps_path, timeout=300):
    """Solve with external MIP solver; return (is_optimal, obj_value) or (False, None)."""
    try:
        r = subprocess.run(
            ["highs", str(mps_path)],
            capture_output=True, text=True, timeout=timeout + 30
        )
        out = r.stdout + r.stderr
        is_opt = "Status            Optimal" in out
        obj = None
        for line in out.splitlines():
            if "Primal bound" in line:
                try:
                    obj = float(line.split()[-1])
                except ValueError:
                    pass
        return is_opt, obj
    except Exception as e:
        print(f"    Solver error: {e}", file=sys.stderr)
        return False, None


# ──────────────────────────────────────────────────────────────────────────────
# Instance definitions
# ──────────────────────────────────────────────────────────────────────────────

def build_instances():
    insts = []

    # 1. n=8, m=2, integer, Cmax — asymmetric machines, moderate setup
    rng = random.Random(42)
    n, m = 8, 2
    proc  = gen_proc(rng, n, m, lo=6, hi=22, asymmetry=True)
    setup = gen_setup(rng, n, m, lo=3, hi=11)
    insts.append(dict(
        name="upms_n8_m2_int_cmax_s42",
        n=n, m=m, proc=proc, setup=setup, weights=None,
        desc="8 jobs, 2 machines, integer times, Cmax, seed=42"
    ))

    # 2. n=8, m=3, integer, WCT — 3 machines, diverse weights
    rng = random.Random(42)
    n, m = 8, 3
    proc  = gen_proc(rng, n, m, lo=5, hi=18, asymmetry=True)
    setup = gen_setup(rng, n, m, lo=2, hi=9)
    w = [rng.randint(1, 5) for _ in range(n)]
    insts.append(dict(
        name="upms_n8_m3_int_wct_s42",
        n=n, m=m, proc=proc, setup=setup, weights=w,
        desc=f"8 jobs, 3 machines, integer, WCT (weights={w}), seed=42"
    ))

    # 3. n=9, m=3, integer, Cmax — harder, different seed
    rng = random.Random(137)
    n, m = 9, 3
    proc  = gen_proc(rng, n, m, lo=7, hi=25, asymmetry=True)
    setup = gen_setup(rng, n, m, lo=3, hi=12)
    insts.append(dict(
        name="upms_n9_m3_int_cmax_s137",
        n=n, m=m, proc=proc, setup=setup, weights=None,
        desc="9 jobs, 3 machines, integer times, Cmax, seed=137"
    ))

    # 4. n=8, m=3, fractional (0.5 grid), Cmax — non-integer optimal
    rng = random.Random(42)
    n, m = 8, 3
    proc  = gen_proc_frac(rng, n, m, lo=4.0, hi=18.0)
    setup = gen_setup_frac(rng, n, m, lo=1.5, hi=8.0)
    insts.append(dict(
        name="upms_n8_m3_frac_cmax_s42",
        n=n, m=m, proc=proc, setup=setup, weights=None,
        desc="8 jobs, 3 machines, fractional proc/setup (0.5 grid), Cmax, seed=42"
    ))

    # 5. n=9, m=2, large setup times (setups ≈ proc times), Cmax
    rng = random.Random(137)
    n, m = 9, 2
    proc  = gen_proc(rng, n, m, lo=5, hi=18, asymmetry=True)
    setup = gen_setup(rng, n, m, lo=8, hi=22)  # large setups
    insts.append(dict(
        name="upms_n9_m2_lgset_cmax_s137",
        n=n, m=m, proc=proc, setup=setup, weights=None,
        desc="9 jobs, 2 machines, large setup times (8-22), Cmax, seed=137"
    ))

    # 6. n=7, m=3, integer, WCT — different seed, certifiable in reasonable time
    rng = random.Random(137)
    n, m = 7, 3
    proc  = gen_proc(rng, n, m, lo=5, hi=18, asymmetry=True)
    setup = gen_setup(rng, n, m, lo=2, hi=9)
    w = [rng.randint(1, 5) for _ in range(n)]
    insts.append(dict(
        name="upms_n7_m3_int_wct_s137",
        n=n, m=m, proc=proc, setup=setup, weights=w,
        desc=f"7 jobs, 3 machines, integer, WCT (weights={w}), seed=137"
    ))

    return insts


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    fixture_dir = Path(__file__).parent / "fixtures"
    fixture_dir.mkdir(exist_ok=True)

    instances = build_instances()
    results = []

    print("=== Generating UPMSP-ST Fixtures ===\n")

    for inst in instances:
        name = inst["name"]
        print(f"[{name}]  {inst['desc']}")

        mps_text = make_upms_mps(
            name, inst["n"], inst["m"],
            inst["proc"], inst["setup"], inst.get("weights")
        )

        mps_path = fixture_dir / f"{name}.mps.gz"
        with gzip.open(mps_path, "wt") as f:
            f.write(mps_text)
        print(f"  Written: {mps_path}  ({mps_path.stat().st_size} bytes)")

        # Verify with external solver
        print("  Verifying with external solver …", end="", flush=True)
        is_opt, obj = verify_with_solver(mps_path, timeout=300)
        if is_opt and obj is not None:
            print(f" optimal={obj:.4f}")
            results.append((name, obj, True))
        else:
            print(f" not proven optimal (obj={obj})")
            results.append((name, obj, False))

    # Summary table for embedding in the C test file
    print("\n=== Summary (paste into C test) ===")
    print("static const UpmsTestCase upms_test_cases[] = {")
    for name, obj, is_opt in results:
        obj_str = f"{obj:.4f}" if obj is not None else "0.0 /* TBD */"
        cert = "/* certified */" if is_opt else "/* not proven */"
        print(f'  {{"{name}", {obj_str}, {cert} 500000, 300}},')
    print("};")

    return 0


if __name__ == "__main__":
    sys.exit(main())
