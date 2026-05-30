---
inclusion: always
---

# MIPster — Product Overview

**MIPster** is an open-source mixed-integer linear programming (MIP) solver written in C++, forked from COIN-OR Cbc (`github.com/h-g-s/mipster`). It is rebranded from Cbc and diverges from upstream in several important ways.

## What was removed from upstream Cbc

- **GLPK/GMPL support** — removed from CoinUtils and Cbc
- **Non-Clp LP solvers** — Cplex, DyLP, Mosek, SoPlex, HiGHS, Symphony removed from Cbc and Cgl
- **GPL dependencies** — readline (interactive mode) and nauty removed
- **Interactive mode** — `mipster` binary no longer enters interactive mode when called without arguments; it prints help and exits
- **Osi project** — merged into the monorepo

## Public API

The primary public API is the **C interface** in `Cbc_C_Interface.cpp/.h`. This is used by downstream projects such as `python-mip`. C++ API compatibility (`CbcMain0`/`CbcMain1`) is not a concern.

## CI & Packaging

GitHub Actions workflows (`.github/workflows/`) build/package for three platforms:

| Workflow | Produces |
|---|---|
| `linux-ci.yml` | Manylinux (x86-64 + aarch64) `.tar.gz`, `.deb`, `.rpm`, `.pkg.tar.zst` |
| `macos-ci.yml` | macOS `.pkg` installer |
| `windows-ci.yml` | Windows NSIS installer |

The library is named **`libmipster`**; headers install to `include/mipster/`.
