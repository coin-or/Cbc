// Copyright (C) 2026, COIN-OR Foundation
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file cbcSolverExample.cpp
    \brief Demonstrates the CbcSolver programmatic API.

    Three usage patterns:
    1. One-call solve: CbcSolver::solve(filename)
    2. Import + configure + solve via run()
    3. Direct parameter control before solving
*/

#include "CbcSolver.hpp"
#include <cstdio>

int main(int argc, const char *argv[])
{
  if (argc < 2) {
    printf("Usage: %s problem.mps [timeLimit]\n", argv[0]);
    return 1;
  }
  const char *filename = argv[1];
  double timeLimit = (argc > 2) ? atof(argv[2]) : 300.0;

  // ── Pattern 1: One-call solve ──────────────────────────────────
  printf("=== Pattern 1: One-call solve ===\n");
  {
    CbcSolver cbc;
    cbc.solve(filename);

    printf("Status: %d  Obj: %g  Nodes: %d  Iters: %d\n",
      cbc.status(), cbc.objectiveValue(),
      cbc.nodeCount(), cbc.iterationCount());
    if (cbc.hasSolution())
      printf("Solution found!\n");
    else
      printf("No solution found.\n");
  }

  // ── Pattern 2: Import + configure + run ────────────────────────
  printf("\n=== Pattern 2: Import + configure + run ===\n");
  {
    CbcSolver cbc;
    cbc.importModel(filename);

    // Set parameters between import and solve
    CbcParameters &params = cbc.parameters();
    params[CbcParam::TIMELIMIT]->setVal(timeLimit);
    params[CbcParam::GOMORYCUTS]->setVal("off");
    params[CbcParam::PROBINGCUTS]->setVal("on");

    // Solve via command queue
    std::deque<std::string> q = {"-solve", "-quit"};
    cbc.run(q);

    printf("Status: %d  Obj: %g  Nodes: %d\n",
      cbc.status(), cbc.objectiveValue(), cbc.nodeCount());
  }

  // ── Pattern 3: Full control with argc/argv ─────────────────────
  printf("\n=== Pattern 3: argc/argv interface ===\n");
  {
    OsiClpSolverInterface solver;
    CbcSolver cbc(solver);
    cbc.initialize();

    // Pass command-line style arguments
    char timeBuf[32];
    snprintf(timeBuf, sizeof(timeBuf), "%.0f", timeLimit);
    const char *args[] = {"cbc", filename, "-sec", timeBuf,
      "-cuts", "off", "-solve", "-quit"};
    cbc.run(8, args);

    printf("Status: %d  Obj: %g  Bound: %g\n",
      cbc.status(), cbc.objectiveValue(), cbc.bestBound());
  }

  return 0;
}
