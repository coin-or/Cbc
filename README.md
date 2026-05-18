# MIPster — Mixed-Integer Programming Solver

**MIPster** is a high-performance open-source solver for Mixed-Integer Linear Programming (MIP/MILP) problems. It is a streamlined fork of [COIN-OR CBC](https://github.com/coin-or/Cbc), redesigned for simplicity, performance, and ease of integration.

---

## Goals

- **Strong default performance** — particularly on problems with set packing, covering, and partitioning constraints, as commonly found in column-generation formulations.
- **Good parallel performance** — opportunistic parallel LP solving via a racing portfolio of LP configurations.
- **Comprehensive C API** — a clean C interface callable from virtually any language, including Python (via CFFI), Julia, and others.
- **Lightweight codebase** — experimental features that do not consistently improve performance on a broad set of instances are removed, not accumulated.
- **Simple build system** — a single `./configure && make && sudo make install`; no multi-project orchestration required.
- **Modern C++** — new code targets C++17.

## Non-Goals

- Supporting multiple LP solver backends (only CLP, already embedded).
- Supporting multiple sparse linear system libraries for the barrier solver (only AMD is kept).

---

## What Changed from Upstream CBC

MIPster started from the CBC 3.x development branch and underwent significant consolidation:

| | Upstream (5 projects) | MIPster |
|---|---|---|
| Projects / build systems | 5 (CoinUtils, Osi, Clp, Cgl, Cbc) | **1** |
| Source lines (C++/headers) | ~572,000 | **~509,000** |
| External dependencies | ~10 optional | **1 optional (AMD)** |
| Solver binary | `cbc` | **`mipster`** |

### Code removed

- **Aboca / AbcSimplex** (~35,000 lines) — experimental parallel simplex code.
- **GLPK / GMPL support** — LP/MIP format parsing via GLPK (GPL-licensed).
- **Multiple Cholesky solvers** — CHOLMOD (GPL), MUMPS, WSMP (proprietary), Pardiso (proprietary), Taucs; only AMD (BSD-licensed) is kept for the barrier method.
- **ClpPdco / ClpLsqr** (~1,900 lines) — interior-point solver for regularized least-squares problems; not used in MIP solving.
- **Standalone `clp` binary** (~4,400 lines) — CLP is now an internal library.
- **Osi abstraction layer** — the solver-independent LP interface is embedded directly; no longer a separate project.
- **All `#if 0` dead-code blocks** — hundreds of commented-out code blocks removed throughout.
- **MSVC-specific workarounds** — Windows builds use MinGW/GCC; MSVC-only code removed.
- **`CoinUtils`, `Clp`, `Cgl` as separate projects** — all merged into a single unified build.

---

## Building

Requirements: a C++17 compiler (GCC or Clang), `make`, and optionally AMD (SuiteSparse) for the barrier LP method.

```sh
./configure
make -j$(nproc)
sudo make install
```

This produces the `mipster` binary and the `libCbc` library.

### Optional: AMD (barrier solver Cholesky)

If SuiteSparse/AMD is installed:
```sh
./configure --with-amd
```

---

## Usage

```sh
# Solve a MPS file
mipster problem.mps -solve

# With a time limit (seconds)
mipster problem.mps -sec 300 -solve

# Write solution to file
mipster problem.mps -solve -solu solution.txt

# Interactive help
mipster --help
```

MIPster reads `.mps`, `.lp`, `.mps.gz`, and `.lp.gz` files directly.

---

## C API

MIPster exposes the same C API as CBC, making it a drop-in replacement:

```c
#include "Cbc_C_Interface.h"

Cbc_Model *model = Cbc_newModel();
Cbc_readMps(model, "problem.mps");
Cbc_solve(model);
printf("Objective: %f\n", Cbc_getObjValue(model));
Cbc_deleteModel(model);
```

---

## License

MIPster is licensed under the [Eclipse Public License 2.0 (EPL-2.0)](LICENSE), the same license as COIN-OR CBC.

---

## Acknowledgements

MIPster is built on top of [COIN-OR CBC](https://github.com/coin-or/Cbc), [CLP](https://github.com/coin-or/Clp), [CGL](https://github.com/coin-or/Cgl), [OSI](https://github.com/coin-or/Osi), and [CoinUtils](https://github.com/coin-or/CoinUtils). We are grateful to the entire COIN-OR community for building and maintaining these foundations.
