// CbcPostprocessRepair.hpp
//
// Repairs the back-substituted integer solution after CglPreProcess::postProcess().
//
// CglPreProcess can transform ≥-rows to ≤-rows (sign-flipping) during
// preprocessing, and its back-substitution uses LP solves to set eliminated
// integer variable values.  Both can produce infeasible integer assignments
// in the original-space solution.  The repair pass uses the original problem
// constraints (from originalSolver) to detect and correct violations through
// three phases:
//
//   Phase 1 – Direction-aware integer rounding for eliminated variables.
//   Phase 2 – Greedy variable-first repair (commits moves with gain > 0).
//   Phase 3 – Tabu search + depth-limited DFS for any remaining violations.
//
// After the repair, saveSolver is updated with the corrected solution and
// integer variable bounds are fixed so that subsequent LP solves preserve
// the repaired values.

#ifndef CbcPostprocessRepair_H
#define CbcPostprocessRepair_H

class OsiSolverInterface;
class CbcModel;
class CglPreProcess;

// Repair infeasible integer variable assignments in saveSolver.
// Modifies saveSolver in-place: sets integer variable bounds and solution
// to a feasible (or best-found) assignment consistent with originalSolver's
// constraint structure.
void CbcRepairPostprocessSolution(
  OsiSolverInterface *saveSolver,
  OsiSolverInterface *originalSolver,
  CbcModel *babModel,
  CglPreProcess &process);

#endif // CbcPostprocessRepair_H
