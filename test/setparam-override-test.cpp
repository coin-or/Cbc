// Parameters set by name through the C interface (Cbc_setParam) must override
// the typed defaults (Cbc_setAllowableFractionGap and friends) that Cbc_solve()
// translates into command-queue strings. Both kinds end up in the same queue,
// and the last setting of a parameter wins, so the typed values have to go in
// front of the ones set by name.
// This file is licensed under the Eclipse Public License (EPL).
//
// Minimize -3x - 2y + 100000z
// s.t.     2x + 3y <= 4
//          x, y binary, z fixed at 1
//
// The optimum is x=1, y=0 (99997). The MIP start x=0, y=1 (99998) is within
// the default ratio gap of 1e-4 of the root bound, so with the default gap the
// search stops at the MIP start, and the optimum is only found when the
// ratioGap set by name takes effect.
//
// Usage: setparam-override-test
// Exit code: 0 = all checks passed, 1 = one or more checks failed.

#include "Cbc_C_Interface.h"

#include <cmath>
#include <cstdio>

namespace {

Cbc_Model *buildModel()
{
  Cbc_Model *model = Cbc_newModel();

  CoinBigIndex start[] = { 0, 1, 2, 2 };
  int rowindex[] = { 0, 0 };
  double value[] = { 2, 3 };
  double collb[] = { 0, 0, 1 };
  double colub[] = { 1, 1, 1 };
  double obj[] = { -3, -2, 100000 };
  double rowlb[] = { -INFINITY };
  double rowub[] = { 4 };
  int mscols[] = { 0, 1 };
  double msvalues[] = { 0, 1 };

  Cbc_loadProblem(model, 3, 1, start, rowindex, value, collb, colub, obj,
    rowlb, rowub);
  Cbc_setInteger(model, 0);
  Cbc_setInteger(model, 1);

  Cbc_setParam(model, "preprocess", "off");
  Cbc_setParam(model, "heur", "off");
  Cbc_setMIPStartI(model, 2, mscols, msvalues);
  return model;
}

int check(bool ok, const char *what)
{
  printf("%s: %s\n", ok ? "ok" : "FAILED", what);
  return ok ? 0 : 1;
}

// Control: with the default ratio gap the search stops at the MIP start.
// If this one fails, the model no longer tells the two cases apart.
int checkDefaultGap()
{
  Cbc_Model *model = buildModel();
  Cbc_solve(model);
  double objValue = Cbc_getObjValue(model);
  printf("default gap: objective %.6f\n", objValue);
  int failures = check(fabs(objValue - 99998.0) < 1e-6,
    "default ratio gap stops at the MIP start (99998)");
  Cbc_deleteModel(model);
  return failures;
}

// ratioGap set by name must replace the typed default gap.
int checkNamedGap()
{
  Cbc_Model *model = buildModel();
  Cbc_setParam(model, "ratioGap", "0");
  Cbc_solve(model);
  double objValue = Cbc_getObjValue(model);
  printf("ratioGap 0 by name: objective %.6f\n", objValue);
  int failures = check(Cbc_isProvenOptimal(model) != 0,
    "ratioGap 0 set by name: proven optimal");
  failures += check(fabs(objValue - 99997.0) < 1e-6,
    "ratioGap 0 set by name reaches the optimum (99997)");
  Cbc_deleteModel(model);
  return failures;
}

} // namespace

int main()
{
  int failures = 0;
  failures += checkDefaultGap();
  failures += checkNamedGap();
  if (failures) {
    printf("setparam-override-test: %d check(s) failed\n", failures);
    return 1;
  }
  printf("setparam-override-test: all checks passed\n");
  return 0;
}
