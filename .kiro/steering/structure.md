---
inclusion: always
---

# MIPster — Repository Structure

## Top-level Layout

```
Cbc/           Main monorepo — BUILD THIS (contains all sources)
  src/         All C++ sources (CoinUtils, Clp, Osi, Cgl, Cbc merged)
  test/        Test suite
  BuildTools/  Autotools m4 macros (symlink)
CoinUtils/     Reference-only (separate upstream project)
Clp/           Reference-only (separate upstream project)
Cgl/           Reference-only (separate upstream project)
Cbc/scripts/configster   All-in-one configure+build+install script
cbc_build.sh       Incremental rebuild script
cbc_clean.sh       Clean script
doc/           Documentation (howtos, debugging guides)
.github/       CI workflows + Copilot instructions (symlink → AGENTS.md)
.kiro/         Kiro steering files
AGENTS.md      Primary AI agent instructions (all tools)
CLAUDE.md      Claude Code entry point (imports AGENTS.md)
```

## Source Organisation within `Cbc/src/`

```
Cbc/src/
  CoinUtils/   Utilities (matrix, messages, clique utilities)
  Clp/         LP solver
  Osi/         Abstract solver interface (OsiSolverInterface, etc.)
  OsiClp/      Clp implementation of Osi
  Cgl/         Cut generators (one subdir per cut family)
  *.cpp/*.hpp  Cbc core (model, solver, heuristics, branching)
```

## Important Non-Source Files

| File | Purpose |
|---|---|
| `Cbc/src/Makefile.am` | Autotools source list — edit when adding new `.cpp`/`.hpp` files |
| `.clang-format` | Code style config (WebKit, 2-space indent) |
| `solvable_instances.tsv` | Instances solvable in <5 min — use for smoke tests |
| `doc/cut-debugging-howto.md` | Guide for debugging invalid cuts |
| `doc/lp_racing_analysis_howto.md` | LP relaxation & racing portfolio analysis workflow |
