#!/usr/bin/env python3
"""Generate Multi-Item Capacitated Lot Sizing with Batch Production (MICLSP) fixtures.

Problem (minimize total setup + holding cost):
  min   sum_{i,t}  f[i]*y[i,t]  +  h[i]*s[i,t]
  s.t.
    inventory balance : s[i,t-1] + B[i]*n[i,t] - s[i,t]  = d[i,t]   ∀i,t
                        (with s[i,-1] = 0)
    machine capacity  : sum_i B[i]*n[i,t]                <= C[t]     ∀t
    setup linking     : n[i,t]                            <= M[i,t]*y[i,t]  ∀i,t
                        (Big-M: M[i,t] = floor(C[t] / B[i]))

  y[i,t] ∈ {0,1}   (binary  – setup indicator)
  n[i,t] ∈ Z+      (integer – number of production batches)
  s[i,t] ∈ R+      (continuous – end-of-period inventory)

Capacity formula (guarantees feasibility for any rho >= 1.0):
  C[t] = ceil( sum_i ceil(d[i,t] / B[i]) * B[i]  *  rho )
  The inner sum is the lot-for-lot capacity (always feasible); rho adds slack.

Instance data generation:
  B[i] ~ U[3,15], d[i,t] ~ U[1,20], f[i] ~ U[100,500], h[i] ~ U[1,10]

Known MIP optima are pre-certified and embedded directly; no re-solve required.

Output: miclsp_fixtures.h
Run:  python3 gen_miclsp_fixtures.py
"""

import sys
import math
import random

# ---------------------------------------------------------------------------
# Instance specs: (I, T, seed, rho, known_optimal)
# Optima certified independently; instances require B&B (not solved at root).
# ---------------------------------------------------------------------------
SPECS = [
    (4, 10,  42, 1.10, 6538),
    (4, 10,  42, 1.15, 6490),
    (4, 10,  42, 1.20, 6395),
    (4, 10,  42, 1.30, 6238),
    (4, 10, 137, 1.10, 8423),
    (4, 10, 137, 1.15, 8302),
    (4, 10, 137, 1.20, 7962),
    (4, 12,  42, 1.10, 8414),
    (4, 12,  42, 1.15, 8316),
    (4, 12,  42, 1.20, 8181),
    (4, 12, 137, 1.15, 8555),
]


def make_instance(I, T, seed, rho):
    rng = random.Random(seed)
    B = [rng.randint(3, 15) for _ in range(I)]
    d = [[rng.randint(1, 20) for _ in range(T)] for _ in range(I)]
    f = [rng.randint(100, 500) for _ in range(I)]
    h = [rng.randint(1, 10) for _ in range(I)]
    # Lot-for-lot capacity * rho: always feasible, tightness controlled by rho.
    C = [int(math.ceil(
             sum(math.ceil(d[i][t] / B[i]) * B[i] for i in range(I)) * rho))
         for t in range(T)]
    return B, d, f, h, C


def arr_c(values, indent=4, per_line=12):
    lines, row = [], []
    for v in values:
        row.append(str(v))
        if len(row) == per_line:
            lines.append(" " * indent + ", ".join(row) + ",")
            row = []
    if row:
        lines.append(" " * indent + ", ".join(row))
    return "\n".join(lines)


def main():
    instances = []
    for I, T, seed, rho, opt in SPECS:
        B, d, f, h, C = make_instance(I, T, seed, rho)
        instances.append((I, T, seed, rho, opt, B, d, f, h, C))

    print("/* Multi-Item Capacitated Lot Sizing with Batch Production (MICLSP)")
    print(" * test fixtures.  All three variable types: binary (setup indicator),")
    print(" * integer (production batches), continuous (inventory).")
    print(" * Big-M links batch count to setup: n[i,t] <= floor(C[t]/B[i]) * y[i,t].")
    print(" * Capacity formula: C[t] = ceil(sum_i ceil(d[i,t]/B[i])*B[i] * rho).")
    print(" *")
    print(" * Do not edit manually; regenerate with test/gen_miclsp_fixtures.py.")
    print(" */")
    print()
    print("#ifndef MICLSP_FIXTURES_H")
    print("#define MICLSP_FIXTURES_H")
    print()
    print("typedef struct {")
    print("  int I;       /* number of items */")
    print("  int T;       /* number of periods */")
    print("  double rho;  /* capacity tightness factor */")
    print("  int opt;     /* certified MIP optimum */")
    print("  const int *B;  /* batch sizes [I] */")
    print("  const int *d;  /* demands [I*T], row-major: d[i*T+t] */")
    print("  const int *f;  /* setup costs [I] */")
    print("  const int *h;  /* holding costs [I] */")
    print("  const int *C;  /* machine capacity [T] */")
    print("} MiclspInstance;")
    print()

    tags = []
    for (I, T, seed, rho, opt, B, d, f, h, C) in instances:
        tag = f"miclsp_i{I:02d}_t{T:02d}_s{seed:03d}_r{int(rho*100):03d}"
        flat_d = [d[i][t] for i in range(I) for t in range(T)]
        print(f"/* {tag}: I={I} T={T} seed={seed} rho={rho:.2f} opt={opt} */")
        print(f"static const int B_{tag}[{I}] = {{")
        print(arr_c(B)); print("};")
        print(f"static const int d_{tag}[{I*T}] = {{")
        print(arr_c(flat_d)); print("};")
        print(f"static const int f_{tag}[{I}] = {{")
        print(arr_c(f)); print("};")
        print(f"static const int h_{tag}[{I}] = {{")
        print(arr_c(h)); print("};")
        print(f"static const int C_{tag}[{T}] = {{")
        print(arr_c(C)); print("};")
        print()
        tags.append((tag, I, T, rho, opt))

    print(f"#define N_MICLSP_INSTANCES {len(tags)}")
    print()
    print("static const MiclspInstance MICLSP_INSTANCES[] = {")
    for (tag, I, T, rho, opt) in tags:
        print(f"  {{ {I}, {T}, {rho:.2f}, {opt},")
        print(f"    B_{tag}, d_{tag}, f_{tag}, h_{tag}, C_{tag} }},")
    print("};")
    print()
    print("#endif /* MICLSP_FIXTURES_H */")


if __name__ == "__main__":
    main()
