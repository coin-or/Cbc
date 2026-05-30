---
inclusion: always
---

# MIPster — Tech Stack & Architecture

## Language & Toolchain

- **C++** (C++14), built with **autotools** (`automake`/`autoconf`/`libtool`)
- Single monorepo (`Cbc/`) embeds CoinUtils, Clp, Osi, Cgl, and Cbc sources under `Cbc/src/`
- **clang-format** (WebKit-based, 2-space indent, no column limit) — see `.clang-format`

## Dependency Order

```
CoinUtils → Clp (includes Osi) → Cgl → Cbc
```

## Project Layers

| Layer | Role |
|---|---|
| **CoinUtils** | Utility types: `CoinPackedMatrix`, `CoinMessage`, `CoinBronKerbosch`, clique utilities |
| **Clp** | LP solver (simplex); includes merged Osi abstract interface (`src/Osi/`, `src/OsiClp/`) |
| **Cgl** | Cut generators: `CglBKClique`, `CglGomory`, `CglOddWheel`, `CglCliqueStrengthening`, etc. |
| **Cbc** | MIP solver: B&B orchestration using Clp as LP relaxation and Cgl cuts |

## Key Cbc Source Files

| File | Purpose |
|---|---|
| `CbcModel.cpp/.hpp` | Core B&B model — central class for the solver |
| `CbcSolver.cpp/.hpp` | Command-line driver and parameter handling |
| `CbcParameters.cpp/.hpp` | Parameter definitions |
| `Cbc_C_Interface.cpp/.h` | Public C API |
| `CbcHeuristic*.cpp` | Primal heuristics (FPump, RINS, RENS, Dive variants) |
| `CbcBranch*.cpp` | Branching objects and strategies |
| `CbcTree.cpp/.hpp` | B&B search tree |
| `CbcInstanceFeatures.cpp/.hpp` | Extracts MIP instance features to CSV |

## Conflict Graph Infrastructure

Recent addition enabling advanced cut generators:
- `CglBKClique` — Bron-Kerbosch clique separator (default, replaces `CglClique`)
- `CglOddWheel` — Odd-wheel cut separator
- `CglCliqueStrengthening` — Preprocessing for set-packing constraints

Controlled via `-cgraph`, `-clqstr`, `-clique`, `-oddwheel` solver parameters.
