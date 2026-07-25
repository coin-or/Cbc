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

// Cheap feasibility check: recomputes every row's activity from saveSolver's
// current solution against originalSolver's row bounds/matrix (the full,
// pre-LP-presolve constraint set) and returns true iff every row is satisfied
// within tolerance. Intended to be called BEFORE CbcRepairPostprocessSolution
// so that its (potentially expensive, and occasionally harmful — see the
// "fiball" postprocess-infeasibility investigation) repair phases are only
// invoked when actually needed.
bool CbcPostprocessSolutionIsFeasible(
  OsiSolverInterface *saveSolver,
  OsiSolverInterface *originalSolver,
  double tolerance = 1.0e-6);

#endif // CbcPostprocessRepair_H
