/* Multi-Echelon Supply Chain Network (MESCN) test.
 *
 * Tests MIP solve correctness for supply chain optimization problems with:
 *   - Binary variables (setup decisions)
 *   - General integer variables (batch counts, truck counts)
 *   - Continuous variables (production, shipments, inventory)
 *   - Flow balance constraints
 *   - Capacity constraints
 *   - Big-M constraints (setup links)
 *
 * Verifies:
 *   1. Objective value matches certified optimal (or correctly identifies infeasible)
 *   2. Solution feasibility (all constraints satisfied)
 *   3. Variable integrality
 *
 * On failure: runs mip_diag_wrong_optimal() to identify which feature is responsible.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mescn_builder.h"
#include "mip_diag.h"
#include "test_utils.h"

/* ── Solution checker ─────────────────────────────────────────────────────
 * Validates a solution vector for semantic correctness.
 * Returns 1 on pass, 0 on any violation.
 */
static int check_mescn_sol(const double *sol, double sol_obj,
                            const MESCNInstance *inst, const char *label)
{
  int F = inst->F, W = inst->W, R = inst->R, P = inst->P, T = inst->T;

  /* Check integrality of binary and integer variables */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        /* u[f,p,t] binary */
        double u_val = sol[uidx(F, W, R, P, T, f, p, t)];
        if (u_val < -1e-6 || u_val > 1.0 + 1e-6 || fabs(u_val - round(u_val)) > 1e-4) {
          printf("    FAIL [%s]: u(%d,%d,%d)=%.6f not binary\n", label, f, p, t, u_val);
          return 0;
        }

        /* n[f,p,t] integer */
        double n_val = sol[nidx(F, W, R, P, T, f, p, t)];
        if (n_val < -1e-6 || fabs(n_val - round(n_val)) > 1e-4) {
          printf("    FAIL [%s]: n(%d,%d,%d)=%.6f not integer\n", label, f, p, t, n_val);
          return 0;
        }
      }
    }
  }

  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        /* v[w,p,t] binary */
        double v_val = sol[vidx(F, W, R, P, T, w, p, t)];
        if (v_val < -1e-6 || v_val > 1.0 + 1e-6 || fabs(v_val - round(v_val)) > 1e-4) {
          printf("    FAIL [%s]: v(%d,%d,%d)=%.6f not binary\n", label, w, p, t, v_val);
          return 0;
        }

        /* m[w,p,t] integer */
        double m_val = sol[midx(F, W, R, P, T, w, p, t)];
        if (m_val < -1e-6 || fabs(m_val - round(m_val)) > 1e-4) {
          printf("    FAIL [%s]: m(%d,%d,%d)=%.6f not integer\n", label, w, p, t, m_val);
          return 0;
        }
      }
    }
  }

  /* Check factory flow balance */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double s_prev = (t > 0) ? sol[sf_idx(F, W, R, P, T, f, p, t - 1)] : 0.0;
        double x_val = sol[xidx(F, P, T, f, p, t)];
        double s_curr = sol[sf_idx(F, W, R, P, T, f, p, t)];
        double shipments = 0.0;
        for (int w = 0; w < W; w++) {
          if (has_fw_arc(inst, f, w))
            shipments += sol[yidx(F, W, P, T, f, w, p, t)];
        }
        double balance = s_prev + x_val - s_curr - shipments;
        if (fabs(balance) > 1e-4) {
          printf("    FAIL [%s]: factory flow balance violated f=%d p=%d t=%d "
                 "(residual=%.6f)\n",
                 label, f, p, t, balance);
          return 0;
        }
      }
    }
  }

  /* Check warehouse flow balance */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double s_prev = (t > 0) ? sol[sw_idx(F, W, R, P, T, w, p, t - 1)] : 0.0;
        double incoming = 0.0;
        for (int f = 0; f < F; f++) {
          if (has_fw_arc(inst, f, w))
            incoming += sol[yidx(F, W, P, T, f, w, p, t)];
        }
        double s_curr = sol[sw_idx(F, W, R, P, T, w, p, t)];
        double outgoing = 0.0;
        for (int r = 0; r < R; r++) {
          if (has_wr_arc(inst, w, r))
            outgoing += sol[zidx(F, W, R, P, T, w, r, p, t)];
        }
        double balance = s_prev + incoming - s_curr - outgoing;
        if (fabs(balance) > 1e-4) {
          printf("    FAIL [%s]: warehouse flow balance violated w=%d p=%d t=%d "
                 "(residual=%.6f)\n",
                 label, w, p, t, balance);
          return 0;
        }
      }
    }
  }

  /* Check demand satisfaction */
  for (int r = 0; r < R; r++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double supply = 0.0;
        for (int w = 0; w < W; w++) {
          if (has_wr_arc(inst, w, r))
            supply += sol[zidx(F, W, R, P, T, w, r, p, t)];
        }
        double demand = inst->demand[(r * P + p) * T + t];
        if (fabs(supply - demand) > 1e-4) {
          printf("    FAIL [%s]: demand not satisfied r=%d p=%d t=%d "
                 "(supply=%.1f demand=%.1f)\n",
                 label, r, p, t, supply, demand);
          return 0;
        }
      }
    }
  }

  /* Check factory capacity */
  for (int f = 0; f < F; f++) {
    for (int t = 0; t < T; t++) {
      double prod = 0.0;
      for (int p = 0; p < P; p++) {
        prod += sol[xidx(F, P, T, f, p, t)];
      }
      if (prod > inst->cap_factory[f] + 1e-4) {
        printf("    FAIL [%s]: factory capacity exceeded f=%d t=%d "
               "(prod=%.1f cap=%d)\n",
               label, f, t, prod, inst->cap_factory[f]);
        return 0;
      }
    }
  }

  /* Check warehouse capacity */
  for (int w = 0; w < W; w++) {
    for (int t = 0; t < T; t++) {
      double inv = 0.0;
      for (int p = 0; p < P; p++) {
        inv += sol[sw_idx(F, W, R, P, T, w, p, t)];
      }
      if (inv > inst->cap_warehouse[w] + 1e-4) {
        printf("    FAIL [%s]: warehouse capacity exceeded w=%d t=%d "
               "(inv=%.1f cap=%d)\n",
               label, w, t, inv, inst->cap_warehouse[w]);
        return 0;
      }
    }
  }

  /* Check batch production: x = batch_size * n */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double x_val = sol[xidx(F, P, T, f, p, t)];
        double n_val = sol[nidx(F, W, R, P, T, f, p, t)];
        double expected = inst->batch_size[p] * n_val;
        if (fabs(x_val - expected) > 1e-4) {
          printf("    FAIL [%s]: batch production violated f=%d p=%d t=%d "
                 "(x=%.1f batch*n=%.1f)\n",
                 label, f, p, t, x_val, expected);
          return 0;
        }
      }
    }
  }

  /* Check setup link: n <= M_f * u */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double n_val = sol[nidx(F, W, R, P, T, f, p, t)];
        double u_val = sol[uidx(F, W, R, P, T, f, p, t)];
        int M_f = inst->cap_factory[f] / inst->batch_size[p];
        if (n_val > M_f * u_val + 1e-4) {
          printf("    FAIL [%s]: setup link violated f=%d p=%d t=%d "
                 "(n=%.1f > M_f*u=%.1f)\n",
                 label, f, p, t, n_val, M_f * u_val);
          return 0;
        }
      }
    }
  }

  /* Check truck shipments: sum_r z = truck_size * m */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double outgoing = 0.0;
        for (int r = 0; r < R; r++) {
          if (has_wr_arc(inst, w, r))
            outgoing += sol[zidx(F, W, R, P, T, w, r, p, t)];
        }
        double m_val = sol[midx(F, W, R, P, T, w, p, t)];
        double expected = inst->truck_size * m_val;
        if (fabs(outgoing - expected) > 1e-4) {
          printf("    FAIL [%s]: truck shipment violated w=%d p=%d t=%d "
                 "(outgoing=%.1f truck*m=%.1f)\n",
                 label, w, p, t, outgoing, expected);
          return 0;
        }
      }
    }
  }

  /* Check truck setup link: m <= M_w * v */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        double m_val = sol[midx(F, W, R, P, T, w, p, t)];
        double v_val = sol[vidx(F, W, R, P, T, w, p, t)];
        int M_w = (inst->cap_warehouse[w] / inst->batch_size[p]) > 0
                      ? (inst->cap_warehouse[w] / inst->batch_size[p])
                      : 1;
        if (m_val > M_w * v_val + 1e-4) {
          printf("    FAIL [%s]: truck setup link violated w=%d p=%d t=%d "
                 "(m=%.1f > M_w*v=%.1f)\n",
                 label, w, p, t, m_val, M_w * v_val);
          return 0;
        }
      }
    }
  }

  /* Check objective value consistency */
  double obj_computed = 0.0;
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        obj_computed += inst->setup_cost_f[f * P + p] * sol[uidx(F, W, R, P, T, f, p, t)];
        obj_computed += inst->prod_cost[f * P + p] * sol[xidx(F, P, T, f, p, t)];
        obj_computed += inst->holding_f[f * P + p] * sol[sf_idx(F, W, R, P, T, f, p, t)];
      }
    }
  }

  for (int f = 0; f < F; f++) {
    for (int w = 0; w < W; w++) {
      if (!has_fw_arc(inst, f, w))
        continue;
      for (int p = 0; p < P; p++) {
        for (int t = 0; t < T; t++) {
          obj_computed += inst->transport_fw[f * W * P + w * P + p] *
                          sol[yidx(F, W, P, T, f, w, p, t)];
        }
      }
    }
  }

  for (int w = 0; w < W; w++) {
    for (int r = 0; r < R; r++) {
      if (!has_wr_arc(inst, w, r))
        continue;
      for (int p = 0; p < P; p++) {
        for (int t = 0; t < T; t++) {
          obj_computed += inst->transport_wr[w * R * P + r * P + p] *
                          sol[zidx(F, W, R, P, T, w, r, p, t)];
        }
      }
    }
  }

  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        obj_computed += inst->setup_cost_w[w * P + p] * sol[vidx(F, W, R, P, T, w, p, t)];
        obj_computed += inst->holding_w[w * P + p] * sol[sw_idx(F, W, R, P, T, w, p, t)];
      }
    }
  }

  if (fabs(sol_obj - obj_computed) > 1.0) {
    printf("    FAIL [%s]: reported obj=%.1f disagrees with computed=%.1f\n",
           label, sol_obj, obj_computed);
    return 0;
  }

  return 1;
}

/* Model builder for diagnostics */
static Cbc_Model *build_for_diag(void *userdata)
{
  return build_mescn((const MESCNInstance *)userdata);
}

int main(void)
{
  int n_failures = 0;

  for (int i = 0; i < N_MESCN_FIXTURES; i++) {
    const MESCNInstance *inst = &MESCN_FIXTURES[i];
    printf("Testing %s ... ", inst->name);
    fflush(stdout);

    Cbc_Model *model = build_mescn(inst);
    Cbc_setLogLevel(model, 0);
    Cbc_setDblParam(model, DBL_PARAM_TIME_LIMIT, 300.0);
    Cbc_solve(model);

    int is_proven_opt = Cbc_isProvenOptimal(model);
    int is_proven_inf = Cbc_isProvenInfeasible(model);
    double obj = Cbc_getObjValue(model);

    /* Check feasibility of best solution and pool (if solution exists) */
    int feasibility_pass = 1;
    if (!is_proven_inf && Cbc_bestSolution(model) != NULL) {
      /* Check best solution */
      const double *bestSol = Cbc_getColSolution(model);
      double maxViolRow = 0.0; int rowIdx = -1;
      double maxViolCol = 0.0; int colIdx = -1;

      if (!Cbc_checkFeasibility(model, bestSol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
        printf("    FAIL: best solution infeasible  "
               "maxViolRow=%.2e (row %d)  maxViolCol=%.2e (col %d)\n",
               maxViolRow, rowIdx, maxViolCol, colIdx);
        feasibility_pass = 0;
      }

      /* Check all pool solutions */
      int nSol = Cbc_numberSavedSolutions(model);
      for (int s = 0; s < nSol; s++) {
        const double *sol = Cbc_savedSolution(model, s);
        double solObj = Cbc_savedSolutionObjective(model, s);
        maxViolRow = 0.0; rowIdx = -1;
        maxViolCol = 0.0; colIdx = -1;

        if (!Cbc_checkFeasibility(model, sol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
          printf("    FAIL: pool solution %d/%d infeasible (obj=%.2f)  "
                 "maxViolRow=%.2e (row %d)  maxViolCol=%.2e (col %d)\n",
                 s+1, nSol, solObj, maxViolRow, rowIdx, maxViolCol, colIdx);
          feasibility_pass = 0;
        }
      }

      if (!feasibility_pass) {
        n_failures++;
      }
    }

    if (inst->certified_optimal == -1) {
      // Expect infeasible
      if (!is_proven_inf) {
        printf("FAIL: claimed feasible, expected infeasible\n");
        mip_diag_wrong_optimal(build_for_diag, (void *)inst, -1, 900);
        n_failures++;
      } else {
        printf("PASS (infeasible)\n");
      }
    } else if (inst->certified_optimal == -999) {
      // Unknown optimal - just check solution is valid if found
      if (Cbc_bestSolution(model) != NULL) {
        const double *sol = Cbc_getColSolution(model);
        if (!check_mescn_sol(sol, obj, inst, inst->name)) {
          n_failures++;
        } else {
          printf("PASS (feasible, obj=%.0f)\n", obj);
        }
      } else {
        printf("PASS (no solution found)\n");
      }
    } else {
      // Expect optimal or at least feasible within time limit
      if (is_proven_opt) {
        // Proven optimal - check objective
        if (fabs(obj - inst->certified_optimal) > 1.0) {
          printf("FAIL: wrong optimal obj=%.0f != certified=%d\n",
                 obj, inst->certified_optimal);
          mip_diag_wrong_optimal(build_for_diag, (void *)inst,
                                 inst->certified_optimal, 900);
          n_failures++;
        } else {
          // Validate solution
          const double *sol = Cbc_getColSolution(model);
          if (!check_mescn_sol(sol, obj, inst, inst->name)) {
            n_failures++;
          } else {
            printf("PASS (obj=%d)\n", inst->certified_optimal);
          }
        }
      } else if (Cbc_bestSolution(model) != NULL) {
        // Found feasible solution but not proven optimal (timeout OK for hard instances)
        const double *sol = Cbc_getColSolution(model);
        if (!check_mescn_sol(sol, obj, inst, inst->name)) {
          n_failures++;
        } else {
          printf("PASS (feasible, obj=%.0f, not proven optimal)\n", obj);
        }
      } else {
        // No solution found and not proven infeasible - this is a failure
        printf("FAIL: no solution found\n");
        n_failures++;
      }
    }

    Cbc_deleteModel(model);
  }

  if (n_failures > 0) {
    printf("\n%d fixture(s) FAILED\n", n_failures);
    return 1;
  }

  printf("\nAll fixtures PASSED\n");
  return 0;
}
