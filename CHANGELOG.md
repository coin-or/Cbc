# Changelog

All notable changes to MIPster are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-06-05

130 commits since v0.1.0. Major themes: distribution (PyPI wheels + native
installers), correctness fixes, new heuristics (FeasibilityJump), LP
auto-tuning, and a much-expanded test/diagnostics toolkit.

### Added — Distribution

- **Python wheels on PyPI** — `pip install mipster` now ships
  self-contained binary wheels for Linux x86_64 / aarch64, macOS x86_64 /
  arm64, and Windows x86_64. AVX2 / Haswell / NEON variants are
  auto-selected at runtime via CPUID. All third-party runtime libraries
  (libgfortran, libquadmath, libstdc++, libwinpthread, zlib, libbz2) are
  bundled — no system dependencies required.
- **Windows NSIS installer** (`mipster-*-setup-windows-x86_64.exe`) —
  installs to Program Files, adds `bin/` to system PATH, registers in
  Add/Remove Programs, and ships Start Menu shortcuts.
- **Linux distro packages** — `.deb`, `.rpm`, and Arch `.pkg.tar.zst`
  packages for both x86_64 and aarch64.
- **macOS installer packages** — signed `.pkg` installers for x86_64 and
  arm64.
- **Self-contained tarballs** with a bundled `THIRD_PARTY_LICENSES.md`
  covering OpenBLAS (BSD-3), GCC runtime libraries (GPLv3 + RLE 3.1),
  zlib, and libbz2.

### Added — Solver features

- **LP auto-tuning** (`-lpMethod=auto`): random-forest based parameter
  selection (`CbcLpParamScorer`) yields a 1.36× geomean speedup on root
  LP solves. Trained on 12-cluster representative MIPLIB instances.
  Replaces the deprecated `-guess` action.
- **Racing LP portfolio** (`-racingLP on`): multiple CLP configurations
  run in parallel; first-to-finish wins. K=2 / K=3 portfolios were
  selected via k-fold cross-validation. Portfolio size is controlled by
  `-threads`.
- **FeasibilityJump primal heuristic**: NNZ-scaled effort budget with
  stall-based termination. Runs at the root and at deep tree nodes
  (every N levels). New parameters: `-heurDecay`, `-treeCutDepth`,
  `-fpumpSkipIfFeasible`. Replaces the legacy `STOP_CUTS_NOW` define.
- **RINS at deep tree nodes**: extends RINS heuristic into the tree
  rather than running only near the root.
- **`-debugCuts` parameter**: per-generator activation of the row-cut
  debugger, rather than a global toggle.

### Added — Diagnostic tools

- **`mipster_diag`**: automated wrong-optimal triage. Runs a focused
  `debugCuts` pass plus an optional 13-config feature-disabling sweep
  (`-scan`) to pinpoint a faulty cut family or bound-propagation step.
  Output names the offending row / column directly. Documented in
  `doc/cut-debugging-howto.md`.
- **`mipster_validate_sol`**: standalone solution feasibility checker —
  validates bounds, integrality, row activities, and the objective value
  against an MPS file.

### Fixed — Correctness (wrong-optimal / unsound cuts)

- `CglFlowCover` wrong-optimal bug; FCNF fixture tests added.
- `CglProbing`: MIP-unsound generated-integer column cuts during
  preprocessing; in-tree row cuts that used node-local bounds.
- `CglKnapsackCover`: two singleton-cut soundness bugs.
- `CglTwoMIR`: cuts that wrongly excluded the optimal integer solution.
- `CglPreProcess::postProcess()`: wrong objective after postsolve when
  the inner LP presolve was active; stale warm-start basis; missing
  oracle-based repair pass; direction-aware objective warning.
- `CbcBoundPropagation`: stale `colType` cache misclassifying binaries.
- `chooseClpBranch`: dropped a fragile assert in favour of a soft skip
  when the LP solution can't be repaired.
- `Cbc_savedSolutionObj` postprocessing bug; corrected integer variable
  values after preprocessing postprocess.
- Heap-use-after-free in `BabHandlerGuard` / `OsiClpSolverInterface`.
- C API random seed silently ignored.
- `fathomMany` integer-rounding correctness after `doubleCheck()`.

### Changed

- **Build system**: new `configster` wrapper at the repo root provides
  a TUI and CLI for the configure → build → install flow, replacing
  hand-rolled invocations. Auto-detects Homebrew binutils on macOS.
- **macOS support**: build process compatibility merged from PR #1.
- **Library naming**: `libCbc` is now `libmipster`; headers install to
  `include/mipster/`; pkg-config name is `mipster`.
- **Heuristics order**: FeasibilityJump runs before FPump; FJ output
  is suppressed in the tree (kept at root only).
- **Tests**: deterministic node-limit-based gating replaces tight
  time limits; solution-pool validation added across all C interface
  tests; `mip_diag` auto-runs on wrong-optimal CI failures.

### Performance

- Faster MPS reading and `OsiFeatures` extraction.
- Heuristic throttling control via `-heurDecay`.
- AVX2 intrinsic include fix in Clp on systems where `config.h` is not
  in scope.

### Documentation

- New `AGENTS.md` development guide.
- Cut debugging how-to: `doc/cut-debugging-howto.md`.
- LP racing analysis how-to: `doc/lp_racing_analysis_howto.md`.
- Integration testing guide: `doc/testing-guide.md`.
- HTML parameter reference auto-deployed to GitHub Pages.
- Doxygen API docs configuration.

### Tests

- New test suites: JSSP, VRP / VRPPD, MESCN, BPC, FCNF fixtures,
  Steiner Tree (directed flow), MICLSP, attfInst1, kantorovich, and a
  43-fixture pure-binary suite.
- All C interface tests now validate the full solution pool, not just
  the best solution.
- Wrong-optimal failures in CI auto-run `mip_diag` for triage.

### Removed / Deprecated

- `-guess` action: deprecated in favour of `-lpMethod=auto`.
- Hard-coded `STOP_CUTS_NOW` define: replaced by `-treeCutDepth`.

## [0.1.0] — 2026-05

Initial public release of MIPster as a streamlined fork of COIN-OR Cbc.
See [the v0.1.0 release on
GitHub](https://github.com/h-g-s/mipster/releases/tag/v0.1.0) for
details.

[0.2.0]: https://github.com/h-g-s/mipster/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/h-g-s/mipster/releases/tag/v0.1.0
