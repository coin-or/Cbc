// CbcPostprocessRepair.hpp
//
// Repairs the back-substituted integer solution after CglPreProcess::postProcess().
//
// CglPreProcess::postProcess() correctly handles all preprocessing transformations
// including variable sign-flips (complementing). However, multi-pass LP
// back-substitution can leave integer variables with slightly fractional values
// due to cascading FP errors, which may cause constraint violations.
//
// This repair pass detects and corrects violations through three phases:
//
//   Phase 1 – Direction-aware integer rounding for fractional integer variables.
//   Phase 2 – Greedy variable-first repair (commits moves with gain > 0).
//   Phase 3 – Tabu search + depth-limited DFS for any remaining violations.
//
// Key invariant: do NOT override postProcess values for retained variables.
// postProcess is the authoritative source for all variable values including
// those for sign-flipped (complemented) retained variables.
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
