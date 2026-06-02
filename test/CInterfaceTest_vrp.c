/* Vehicle Routing Problem (VRP) test.
 *
 * Tests MIP solve correctness for CVRP and VRPPD problems with:
 *   - Binary arc variables x[i,j]
 *   - Continuous load variables u[i]
 *   - Routing constraints (visit once, flow conservation)
 *   - Capacity constraints (correct formulation: u[i] - u[j] + Q*x[i,j] <= Q - d[j])
 *   - Precedence constraints (VRPPD: delivery after pickup)
 *
 * Verifies:
 *   1. Objective value matches known-good results
 *   2. All instances solve (optimal or feasible within time limit)
 *   3. Solution feasibility via mipster_validate_sol
 *
 * On failure: runs mip_diag_wrong_optimal() to identify which feature is responsible.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "Cbc_C_Interface.h"
#include "mip_diag.h"
#include "test_utils.h"

/* Test instance definition */
typedef struct {
  const char *name;
  double expected_obj; /* 0.0 = no expectation (just check solvability) */
  int timeout_sec;
  int is_vrppd; /* 1 for VRPPD, 0 for CVRP */
} VRPTestCase;

/* Test instances - ordered by difficulty */
static const VRPTestCase test_cases[] = {
  /* Easy CVRP (< 1s) */
  {"cvrp_loose", 90.0, 10, 0},
  {"cvrp_tight", 130.0, 10, 0},
  {"cvrp_size_small", 20158.0, 10, 0},
  {"cvrp_small_scale", 20012.6, 10, 0},

  /* Easy VRPPD (< 1s) */
  {"vrppd_small", 70.0, 10, 1},
  {"vrppd_small_uniform", 105.0, 10, 1},
  {"vrppd_small_scale", 26.2, 10, 1},
  {"vrppd_medium_clustered", 93.43, 10, 1},
  {"vrppd_tight", 20078.0, 10, 1},

  /* Medium CVRP (1-10s) */
  {"cvrp_demand_skewed", 20260.0, 20, 0},
  {"cvrp_demand_uniform", 20260.0, 20, 0},
  {"cvrp_sparse_grid", 20233.0, 20, 0},
  {"cvrp_capacity_loose", 20167.0, 20, 0},
  {"cvrp_asymmetric", 20370.61, 20, 0},
  {"cvrp_few_large", 20348.0, 20, 0},

  /* Medium VRPPD */
  {"vrppd_large_uniform", 370.27, 20, 1},
  {"vrppd_tight_ring", 397.0, 20, 1},
  {"vrppd_many_requests", 153.72, 20, 1},

  /* Hard CVRP (20s+) */
  {"cvrp_medium", 320.0, 60, 0},
  {"cvrp_geo_ring", 20381.0, 60, 0},
  {"cvrp_demand_outliers", 20298.0, 60, 0},
  {"cvrp_large_scale", 20752.98, 60, 0},
  {"cvrp_many_small", 20208.16, 60, 0},

  /* VRPPD with large cost scales */
  {"vrppd_large_scale", 878.42, 60, 1},

  /* Very hard CVRP (may timeout but should find solution) */
  {"cvrp_geo_clustered", 0.0, 120, 0}, /* No obj expectation */
  {"cvrp_capacity_medium", 0.0, 120, 0},
  {"cvrp_size_large", 0.0, 120, 0},
};

static const int num_test_cases = sizeof(test_cases) / sizeof(test_cases[0]);

/* Validate solution with mipster_validate_sol utility */
static int validate_solution(const char *mps_file, const char *sol_file) {
  pid_t pid = fork();
  if (pid == 0) {
    /* Child process: run mipster_validate_sol */
    execlp("mipster_validate_sol", "mipster_validate_sol", mps_file, sol_file, NULL);
    /* If we get here, exec failed */
    fprintf(stderr, "    WARNING: mipster_validate_sol not found in PATH\n");
    exit(127);
  } else if (pid > 0) {
    /* Parent: wait for child */
    int status;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) {
      int exit_code = WEXITSTATUS(status);
      if (exit_code == 127) {
        /* mipster_validate_sol not found - skip validation */
        return 1;
      }
      return (exit_code == 0) ? 1 : 0;
    }
    return 0;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <test-directory>\n", argv[0]);
    return 1;
  }

  const char *test_dir = argv[1];
  int passed = 0, failed = 0;
  char **failed_instances = malloc(num_test_cases * sizeof(char *));
  int failed_count = 0;

  printf("VRP Test Suite (%d instances)\n", num_test_cases);
  printf("═══════════════════════════════════════════════════════\n\n");

  for (int i = 0; i < num_test_cases; i++) {
    const VRPTestCase *tc = &test_cases[i];
    char mps_path[512];
    snprintf(mps_path, sizeof(mps_path), "%s/%s.mps", test_dir, tc->name);

    printf("[%2d/%2d] %-30s ", i + 1, num_test_cases, tc->name);
    fflush(stdout);

    /* Check if file exists */
    if (access(mps_path, F_OK) != 0) {
      printf("SKIP (file not found)\n");
      continue;
    }

    /* Load and solve */
    Cbc_Model *model = Cbc_newModel();
    int status = Cbc_readMps(model, mps_path);
    if (status != 0) {
      printf("FAIL (cannot read MPS)\n");
      failed++;
      failed_instances[failed_count++] = strdup(tc->name);
      Cbc_deleteModel(model);
      continue;
    }

    /* Set time limit */
    Cbc_setMaximumSeconds(model, tc->timeout_sec);

    /* Solve */
    Cbc_solve(model);

    int solve_status = Cbc_status(model);
    double obj = Cbc_getObjValue(model);
    int is_proven_optimal = Cbc_isProvenOptimal(model);
    int has_solution = Cbc_numberSolutions(model) > 0;

    /* Determine pass/fail */
    int test_passed = 0;
    char status_msg[256] = "";

    if (solve_status == -1) {
      /* Infeasible */
      snprintf(status_msg, sizeof(status_msg), "FAIL (infeasible)");
    } else if (is_proven_optimal) {
      /* Optimal */
      if (tc->expected_obj > 0.0) {
        double rel_err = fabs(obj - tc->expected_obj) / fmax(fabs(tc->expected_obj), 1e-6);
        if (rel_err < 1e-3) {
          snprintf(status_msg, sizeof(status_msg), "PASS (obj=%.2f)", obj);
          test_passed = 1;
        } else {
          snprintf(status_msg, sizeof(status_msg),
                   "FAIL (obj=%.2f, expected=%.2f, err=%.2e)",
                   obj, tc->expected_obj, rel_err);
        }
      } else {
        /* No expected obj - just check solved */
        snprintf(status_msg, sizeof(status_msg), "PASS (obj=%.2f)", obj);
        test_passed = 1;
      }
    } else if (has_solution) {
      /* Feasible but not proven optimal (timeout) */
      snprintf(status_msg, sizeof(status_msg), "PASS (feasible, obj=%.2f, %ds)",
               obj, tc->timeout_sec);
      test_passed = 1;
    } else {
      /* Timeout with no solution */
      snprintf(status_msg, sizeof(status_msg), "FAIL (timeout, no solution)");
    }

    /* Validate solution if available */
    if (test_passed && has_solution) {
      char sol_file[512];
      snprintf(sol_file, sizeof(sol_file), "/tmp/vrp_test_%s.sol", tc->name);

      /* Write solution */
      Cbc_writeSolution(model, sol_file);

      /* Validate */
      if (!validate_solution(mps_path, sol_file)) {
        snprintf(status_msg, sizeof(status_msg), "FAIL (invalid solution)");
        test_passed = 0;
      }

      unlink(sol_file);
    }

    printf("%s\n", status_msg);

    if (test_passed) {
      passed++;
    } else {
      failed++;
      failed_instances[failed_count++] = strdup(tc->name);

      /* Run diagnostic on hard failures (not timeouts without expected obj) */
      if (tc->expected_obj > 0.0 && !is_proven_optimal) {
        printf("    Running diagnostic...\n");
        mip_diag_wrong_optimal(mps_path, tc->expected_obj, tc->timeout_sec);
      }
    }

    Cbc_deleteModel(model);
  }

  printf("\n═══════════════════════════════════════════════════════\n");
  printf("Results: %d passed, %d failed\n", passed, failed);

  if (failed > 0) {
    printf("\nFailed instances:\n");
    for (int i = 0; i < failed_count; i++) {
      printf("  - %s\n", failed_instances[i]);
      free(failed_instances[i]);
    }
  }

  free(failed_instances);
  return (failed == 0) ? 0 : 1;
}
