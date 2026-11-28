/*
 * row-activity-test — targeted tests for Cbc_getRowActivity()/
 * Cbc_getRowSlack() after Cbc_solve(), covering both the pure-LP dispatch
 * (no integers/SOS: Cbc_solve() forwards to Cbc_solveLinearProgram()) and a
 * genuine MIP solve (branch and bound through a CbcModel).
 *
 * This is the regression test for the "rSlk/rActv null/stale pointer after
 * MIP solve" bug: Cbc_getMIPOptimizationResults() used to hand
 * Cbc_updateSlack() cbcModel.getRowActivity() -- which is
 * cbcModel's *internal working solver*'s row activity, i.e. whatever LP was
 * last resolved during the search. That is not necessarily the node that
 * produced the incumbent (it can be a fractional sibling node explored
 * after the best solution was found, or the state the search happened to
 * stop in), so model->rActv/model->rSlk could silently end up describing a
 * *different point* than model->x (the incumbent returned by
 * Cbc_getColSolution()) -- wrong values rather than a crash, and therefore
 * not something an objective-value check would ever catch.
 *
 * Every check below recomputes the expected row activity independently, as
 * A * Cbc_getColSolution() using the row's own coefficients (not by trusting
 * any other Cbc-reported quantity), and confirms Cbc_getRowActivity()/
 * Cbc_getRowSlack() agree with that -- for the *reported* solution, not for
 * whatever the last LP relaxation happened to be. The MIP fixture is chosen
 * so the two disagree on a real build with the bug present: its incumbent
 * fixes z=1 with x=9.5 (row activity x-10z = -0.5), while a node explored
 * later in the search (e.g. proving z=1,x=9.5 optimal after visiting a
 * fractional sibling) can leave the internal solver's LP at a different
 * vertex (e.g. z=0.95 => activity 0), so a stale-activity implementation
 * reports 0 instead of -0.5.
 *
 * Usage: row-activity-test
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "Cbc_C_Interface.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool cond, const std::string &msg)
{
  if (!cond) {
    ++g_failures;
    printf("  FAIL: %s\n", msg.c_str());
  } else {
    printf("  ok:   %s\n", msg.c_str());
  }
}

void checkClose(double actual, double expected, double tol, const std::string &msg)
{
  double delta = fabs(actual - expected);
  if (delta > tol) {
    ++g_failures;
    printf("  FAIL: %s (actual=%.10g expected=%.10g delta=%.3e > tol=%.3e)\n",
      msg.c_str(), actual, expected, delta, tol);
  } else {
    printf("  ok:   %s (actual=%.10g expected=%.10g)\n",
      msg.c_str(), actual, expected);
  }
}

/* A row stated as (indices, coefficients, sense, rhs), matching Cbc_addRow's
 * own argument order, so that the expected activity for a given solution can
 * be recomputed by hand from the same numbers used to build the row. */
struct Row {
  std::string name;
  std::vector< int > idx;
  std::vector< double > coef;
  char sense;
  double rhs;
};

/* A(row, sol): the row's own activity at an arbitrary solution vector,
 * computed directly from (idx, coef) -- independent of anything Cbc reports,
 * so it can serve as ground truth for both the LP and the MIP fixture. */
double rowActivityAt(const Row &row, const double *sol)
{
  double a = 0.0;
  for (size_t k = 0; k < row.idx.size(); ++k)
    a += row.coef[k] * sol[row.idx[k]];
  return a;
}

/* Expected slack, following the same L/G/E convention as Cbc_updateSlack(). */
double expectedSlack(const Row &row, double activity)
{
  switch (row.sense) {
    case 'L':
      return row.rhs - activity;
    case 'G':
      return activity - row.rhs;
    default:
      return fabs(activity - row.rhs);
  }
}

/* Checks Cbc_getRowActivity()/Cbc_getRowSlack() against the ground-truth
 * activity/slack of every row at Cbc_getColSolution(), after model has
 * already been solved. */
void checkRowActivityAndSlack(Cbc_Model *model, const std::vector< Row > &rows,
  const std::string &label)
{
  const double *sol = Cbc_getColSolution(model);
  check(sol != NULL, label + ": Cbc_getColSolution() is non-NULL");
  if (!sol)
    return;

  const double *act = Cbc_getRowActivity(model);
  check(act != NULL, label + ": Cbc_getRowActivity() is non-NULL");
  const double *slk = Cbc_getRowSlack(model);
  check(slk != NULL, label + ": Cbc_getRowSlack() is non-NULL");
  if (!act || !slk)
    return;

  for (size_t i = 0; i < rows.size(); ++i) {
    double expectedActivity = rowActivityAt(rows[i], sol);
    checkClose(act[i], expectedActivity, 1e-6,
      label + ": " + rows[i].name + " activity matches A*x at reported solution");
    checkClose(slk[i], expectedSlack(rows[i], expectedActivity), 1e-6,
      label + ": " + rows[i].name + " slack matches A*x at reported solution");
  }
}

/* ------------------------------------------------------------------ tests */

/* Cbc_solve() on a model with no integers/SOS dispatches straight to
 * Cbc_solveLinearProgram(), skipping the CbcModel/branch-and-bound machinery
 * entirely. Confirms that path's row activity/slack are reported for the LP
 * solution actually returned. */
void testLinearProgram()
{
  printf("\n-- pure LP solve (Cbc_solve() -> Cbc_solveLinearProgram()) --\n");

  Cbc_Model *model = Cbc_newModel();
  Cbc_addCol(model, "x", 0.0, 1e30, 0.0, 0, 0, NULL, NULL);
  Cbc_addCol(model, "y", 0.0, 1e30, 0.0, 0, 0, NULL, NULL);
  Cbc_setObjCoeff(model, 0, 3.0);
  Cbc_setObjCoeff(model, 1, 2.0);

  std::vector< Row > rows;
  {
    Row r{"c1", {0, 1}, {1.0, 1.0}, 'L', 10.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }
  {
    Row r{"c2", {0, 1}, {1.0, -1.0}, 'G', -2.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }
  {
    Row r{"c3", {0}, {1.0}, 'L', 6.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }

  Cbc_setObjSense(model, -1); // maximize
  Cbc_setLogLevel(model, 0);
  Cbc_solve(model);

  check(Cbc_isProvenOptimal(model) != 0, "LP: proven optimal");
  checkRowActivityAndSlack(model, rows, "LP");

  Cbc_deleteModel(model);
}

/* Genuine MIP solve, exercising Cbc_getMIPOptimizationResults()'s row
 * activity/slack path (Cbc_updateSlack()+model->rActv/model->rSlk). x is
 * continuous, y is a general integer, z is binary; c1 (x - 10z <= 0) is the
 * row that discriminates a stale internal-solver activity from the
 * incumbent's own: the optimum fixes z=1 and x=9.5, so c1's activity must
 * come back as -0.5 (slack 0.5), not 0 (which is what a fractional z close
 * to 0.95 at some other explored node would give). */
void testMip()
{
  printf("\n-- genuine MIP solve (Cbc_solve() -> branch and bound) --\n");

  Cbc_Model *model = Cbc_newModel();
  Cbc_addCol(model, "x", 0.0, 1e30, 0.0, 0, 0, NULL, NULL);
  Cbc_addCol(model, "y", 0.0, 1e30, 0.0, 1, 0, NULL, NULL);
  Cbc_addCol(model, "z", 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
  Cbc_setObjCoeff(model, 0, 4.0);
  Cbc_setObjCoeff(model, 1, 1.0);
  Cbc_setObjCoeff(model, 2, -1.0);

  std::vector< Row > rows;
  {
    Row r{"c1", {0, 2}, {1.0, -10.0}, 'L', 0.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }
  {
    Row r{"c2", {0}, {1.0}, 'L', 9.5};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }
  {
    Row r{"c3", {0, 1}, {1.0, 1.0}, 'L', 20.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }

  Cbc_setObjSense(model, -1); // maximize
  Cbc_setLogLevel(model, 0);
  Cbc_solve(model);

  check(Cbc_isProvenOptimal(model) != 0, "MIP: proven optimal");

  const double *sol = Cbc_getColSolution(model);
  if (sol) {
    checkClose(sol[0], 9.5, 1e-6, "MIP: x at expected optimum");
    checkClose(sol[1], 10.0, 1e-6, "MIP: y at expected optimum");
    checkClose(sol[2], 1.0, 1e-6, "MIP: z at expected optimum");
  }

  checkRowActivityAndSlack(model, rows, "MIP");

  Cbc_deleteModel(model);
}

/* Same MIP fixture, but with the *extra* solution pool explicitly disabled
 * (INT_PARAM_MAX_SAVED_SOLS = 0). CbcModel::numberSavedSolutions() always
 * counts the incumbent itself on top of that extra pool (returning
 * numberSavedSolutions_ + 1 whenever a best solution exists), so it comes
 * back as 1 rather than 0 here -- Cbc_getMIPOptimizationResults() therefore
 * still takes the path that populates model->x/model->rActv/model->rSlk,
 * and this exercises that with the smallest possible pool size, rather than
 * with python-mip's default of 10. */
void testMipMinimalSolutionPool()
{
  printf("\n-- genuine MIP solve, extra solution pool disabled --\n");

  Cbc_Model *model = Cbc_newModel();
  Cbc_addCol(model, "x", 0.0, 1e30, 0.0, 0, 0, NULL, NULL);
  Cbc_addCol(model, "y", 0.0, 1e30, 0.0, 1, 0, NULL, NULL);
  Cbc_addCol(model, "z", 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
  Cbc_setObjCoeff(model, 0, 4.0);
  Cbc_setObjCoeff(model, 1, 1.0);
  Cbc_setObjCoeff(model, 2, -1.0);

  std::vector< Row > rows;
  {
    Row r{"c1", {0, 2}, {1.0, -10.0}, 'L', 0.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }
  {
    Row r{"c2", {0}, {1.0}, 'L', 9.5};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }
  {
    Row r{"c3", {0, 1}, {1.0, 1.0}, 'L', 20.0};
    Cbc_addRow(model, r.name.c_str(), (int)r.idx.size(), r.idx.data(), r.coef.data(), r.sense, r.rhs);
    rows.push_back(r);
  }

  Cbc_setIntParam(model, INT_PARAM_MAX_SAVED_SOLS, 0);
  Cbc_setObjSense(model, -1); // maximize
  Cbc_setLogLevel(model, 0);
  Cbc_solve(model);

  check(Cbc_isProvenOptimal(model) != 0, "MIP (minimal pool): proven optimal");
  check(Cbc_numberSavedSolutions(model) == 1,
    "MIP (minimal pool): numberSavedSolutions() is 1 (just the incumbent)");
  checkRowActivityAndSlack(model, rows, "MIP (minimal pool)");

  Cbc_deleteModel(model);
}

} // end anonymous namespace

int main(int argc, char **argv)
{
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      printf("row-activity-test — regression test for Cbc_getRowActivity()/\n"
             "Cbc_getRowSlack() after Cbc_solve(), for both the pure-LP dispatch\n"
             "path and a genuine MIP (branch and bound) solve.\n\n"
             "Usage: row-activity-test\n\n"
             "Takes no options and needs no data files: every fixture and its\n"
             "expected row activity/slack is built and recomputed in-process.\n\n"
             "Exit code: 0 = all checks passed, 1 = one or more failed.\n");
      return 0;
    }
  }

  printf("=== row-activity-test: Cbc_getRowActivity()/Cbc_getRowSlack() after Cbc_solve() ===\n");

  testLinearProgram();
  testMip();
  testMipMinimalSolutionPool();

  printf("\n");
  if (g_failures) {
    printf("RESULT: %d check(s) FAILED\n", g_failures);
    return 1;
  }
  printf("RESULT: all checks passed\n");
  return 0;
}
